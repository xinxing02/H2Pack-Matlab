/*
 * [R, p, r] = QRpartial(A, type, par)
 *
 * Partial pivoting QR decompostion simplified specifically for 
 * interpolative decomposition where the orthogonal matrix is not needed. 
 * 
 * The partial pivoting QR decomposition is of form
 * 
 *      A(:, p) = Q * [R11, R12; 0, R22]
 *
 * where R11 is an upper-triangular matrix and R22 is a dense matrix. 
 * 
 * Input: 
 *      A, target matrix
 *   type, 'rank' or 'tol'
 *    par, parameter for the corresponding rank decision 'type'
 *
 * The dimension of R11 above, is decided by either of the following way 
 *    type = 'rank', par = integer k
 *    type = 'tol',  par = absolute threshold tol.
 *
 * Output: 
 *      R, full matrix [R11, R12; 0, R22]
 *      p, permutation indices shown above
 *      r, dimension of upper-triangular matrix R11. 
 * 
 * 
 * compile command:
 *  mex -O -largeArrayDims QRpartial.c -lmwblas -lmwlapack
 */

#include "mex.h"
#include "matrix.h"
#include "lapack.h"
#include "blas.h"
#include "string.h"
#include "math.h"

void swap_double(double *x, double *y, size_t N)
{
    double tmp; 
    for(size_t i = 0; i < N; i++){
        tmp = *(x+i);
        *(x+i) = *(y+i);
        *(y+i) = tmp;
    }
}

void QRpartial_without_Q(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //  temporary variable for fortran usage 
    ptrdiff_t la, lb, lc, lda, ldb, ldc, ica, icb, icc;
    double scala, scalb, scalc; 
    char trans;

    //  input variable initialization 
    double *A;
    char *type;
    double para; 
    mwSignedIndex m, n, min_mn, n_blk;    
    mwSignedIndex i, j, k;

    //  read the matrix and its dimension
    A = mxGetPr(prhs[0]);   // m * n matrix
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    min_mn = m < n ? m : n;
    
    //  read the stopping criteria and dimensions
    type = mxArrayToString(prhs[1]);
    para = mxGetScalar(prhs[2]);
    if (strcmp(type, "rank") == 0)
        para = para < min_mn ? para : min_mn;

    //  absolute tolerance for each column
    double coltol = sqrt(m) * 1e-15;
    if (strcmp(type, "tol") == 0)
        coltol = (coltol > para) ? coltol : para;

    //  output initialization (R, p)
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, n, mxREAL);

    //  calculation variable initialization
    double *R = mxGetPr(plhs[0]);   //  upper triangle matrix
    double *p = mxGetPr(plhs[1]);   //  permutation vector 
    double *colnorm = (double*)mxMalloc(n * sizeof(double));
    double *h  = (double*)mxMalloc(m * sizeof(double));    // Householder vectors
    double *mv = (double*)mxMalloc(n * sizeof(double));    // h' * A  
    int rank = -1;   
    
    //  input of mex function is read-only variable
    memcpy(R, A, m*n*sizeof(double));
    for (i=0; i<n; i++)
        p[i] = i+1;        

    //  norm of each column
    la = m; ica = 1;
    for (i=0; i<n; i++)
        colnorm[i] = dnrm2(&la, R+i*m, &ica);


    //  iteration of the QR decomposition 
    mwSignedIndex max_iter = min_mn;
    for (i=0; i<max_iter; i++)
    {
        //  Find the pivot column
        mwSignedIndex pivot = i;
        double val = colnorm[i];
        for(j = i; j < n; j++)
            if (colnorm[j] > val){
                val = colnorm[j];
                pivot = j;
            }

        //mexPrintf("%d, %e, %e\n", i, val, para);
        
        //  Check the termination criteria
        if (strcmp(type, "rank") == 0 && (i >= para || val < coltol)){
            rank = i;
            break;
        }
        if (strcmp(type, "tol") == 0 && (val < para || val < coltol)){
            rank = i;   // possibility of stopping before the very first iteration, thus rank equals 0
            break;
        }

        //  Swapping
        if (i != pivot){
            swap_double(colnorm + i, colnorm + pivot, 1);
            swap_double(p + i,  p + pivot, 1);
            swap_double(R + m*i, R + m*pivot, m);
        }

        //  Householder Orthogonalization
        //  vector v for Householder transformation
        memcpy(h, R+m*i+i, (m-i)*sizeof(double));
        int sign = (*h >= 0) - (*h < 0);
        double v_nrm1, v_nrm2;
        la = m-i; ica = 1;  
        v_nrm1 = dnrm2(&la, h, &ica); 
        h[0] = h[0] - (-sign * v_nrm1);
        v_nrm2 = dnrm2(&la, h, &ica);

        //  normalization of h
        la = m-i; scala = 1.0/v_nrm2; ica = 1;
        dscal(&la, &scala, h, &ica);    
     
        //  elimination of the ith sub-column
        *(R+m*i+i) = -sign * v_nrm1;
        memset(R+m*i+(i+1), 0, (m-i-1) * sizeof(double));   


        //  orthogonalization of columns right to the ith column
        trans = 'T';
        la = m-i; lb = n-i-1;
        scala = 1.0; scalb = 0.0; 
        lda = m; icb = 1; icc = 1;
        dgemv(&trans, &la, &lb, &scala, R + m*(i+1) + i, &lda, h, &icb, &scalb, mv, &icc);

        scala = -2.0; 
        dger(&la, &lb, &scala, h, &icb, mv, &icc, R + m*(i+1) + i, &lda);

        //  update the column norm
        for (j=i+1; j<n; j++)
        {
            //  skip small columns
            if (colnorm[j] < coltol)
            {
                colnorm[j] = 0;
                continue; 
            }

            //  update column norm
            double tmp = colnorm[j]*colnorm[j] - (*(R+m*j+i)) * (*(R+m*j+i));
            if (tmp <= 1e-10)
            {
                la = m-i-1; ica = 1;
                colnorm[j] = dnrm2(&la, R+m*j+i+1, &ica);
            }
            else
                colnorm[j] = sqrt(tmp);
        }


        //  check full rank
        if (i == max_iter-1)
        {
            rank = max_iter;
            break;
        }
    }

    if (rank < 0){
    	mexErrMsgTxt("Rank not modified\n");
    	return ;
    }

    //  mexPrintf("%f, %f\n", t1, t2);

    //  output rank
    plhs[2] = mxCreateDoubleScalar(rank);

    //  free space
    mxFree(colnorm);
    mxFree(h);
    mxFree(mv);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* check for proper number of arguments */
    if (nrhs != 3) {
        mexErrMsgTxt("QRpartial() requires 3 arguments.");
    }
    if (!mxIsDouble(prhs[0]) || mxIsSparse(prhs[0])) {
        mexErrMsgTxt( "Input matrix must be a dense matrix." );
    }
    if (!mxIsChar(prhs[1])){
        mexErrMsgTxt( "Stopping criteria must be a string as \"tol\" or \"rank\"." );
    }
    if (!mxIsScalar(prhs[2])){
        mexErrMsgTxt( "Stopping criteria parameter must be a scalar.");
    }

    if (nlhs == 3) {
        QRpartial_without_Q(nlhs, plhs, nrhs, prhs);
    }
    else {
        mexErrMsgTxt( "QR partial with Q component is not supported yet." );
    }
}
