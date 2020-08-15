/*
 * Batched QR decompostion with pivoting
 *
 * [R,p] = QRdim(A, dim, type, par)
 * [Q,R,p] = QRdim(A, dim, type, par)
 *
 * compile command:
 * mex -O QRdim.c libmwblas.lib libmwlapack.lib
 * mex -O -largeArrayDims QRdim_c.c -lmwblas -lmwlapack
 *
 */



#include "mex.h"
#include "matrix.h"
#include "lapack.h"
#include "blas.h"
#include "string.h"
#include "math.h"


inline void swap_int(int *x, int *y, size_t N)
{
    int tmp; 
    for(size_t i = 0; i < N; i++){
        tmp = *(x+i);
        *(x+i) = *(y+i);
        *(y+i) = tmp;
    }
}

inline void swap_double(double *x, double *y, size_t N)
{
    double tmp; 
    for(size_t i = 0; i < N; i++){
        tmp = *(x+i);
        *(x+i) = *(y+i);
        *(y+i) = tmp;
    }
}


void QRdim_without_Q(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //  temporary variable for fortran usage 
    ptrdiff_t la, lb, lc, lda, ldb, ldc, ica, icb, icc;
    double scala, scalb, scalc; 
    char trans;

    //  input variable initialization 
    double *A;
    int dim;
    char *type;
    double para; 
    mwSignedIndex m, n, min_mn, n_blk;    
    mwSignedIndex i, j, k, l;

    //  read the matrix and its dimension
    A = mxGetPr(prhs[0]);
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    min_mn = m < n ? m : n;
    
    //  read the stopping criteria and dimensions
    dim = mxGetScalar(prhs[1]);
    type = mxArrayToString(prhs[2]);
    para = mxGetScalar(prhs[3]);
    n_blk = n / dim;
    if (strcmp(type, "rank") == 0)
        para = para < min_mn ? para : min_mn;

    //  absolute tolerance for each column
    double coltol = sqrt(m) * 1e-15;
    if (strcmp(type, "tol") == 0)
    {
        coltol = (coltol > para) ? coltol : para;
        para   = sqrt((double)dim) * para;
        coltol = sqrt((double)dim) * coltol;
    }
    
    //  output initialization (R, p)
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, n, mxREAL);
    
    //  calculation variable initialization
    double *R = mxGetPr(plhs[0]);   //  upper triangle matrix
    double *p = mxGetPr(plhs[1]);   //  permutation vector 
    double *blknorm = (double*)mxMalloc(n_blk * sizeof(double));
    double *h  = (double*)mxMalloc(m * sizeof(double));   // Householder vectors
    double *mv = (double*)mxMalloc(n * sizeof(double));   // A * h 
    int rank = -1;   
    
    //  input of mex function is read-only variable
    memcpy(R, A, m*n*sizeof(double));
    for (i=0; i<n; i++)
        p[i] = i+1;    
    
    //  norm of each column block (dim columns)
    la = dim * m; ica = 1;
    for(i = 0; i < n_blk; i++)
        blknorm[i] = dnrm2(&la, R + i*m*dim, &ica);
    

    //  iteration of the QR decomposition 
    mwSignedIndex max_iter = min_mn / dim;
    for(i = 0; i < max_iter; i++)
    {
        //  Find the pivot column block
        mwSignedIndex pivot = i;
        double val = blknorm[i];
        for(j = i; j < n_blk; j++)
            if (blknorm[j] > val){
                val = blknorm[j];
                pivot = j;
            }

        //mexPrintf("%d, %e, %e\n", i, val, para);
        
        //  Check the termination criteria
        if (strcmp(type, "rank") == 0 && (i*dim >= para || val < coltol))
        {
            rank = i * dim;
            break;
        }
        if (strcmp(type, "tol") == 0 && (val < para || val < coltol))
        {
            rank = i * dim;     // possibility of stopping before the very first iteration, thus rank equals 0
            break;
        }

        //  Swapping
        if (i != pivot)
        {
            swap_double(blknorm + i, blknorm + pivot, 1);
            swap_double(p + i*dim,  p + pivot*dim, dim);
            swap_double(R + m*i*dim, R + m*pivot*dim, dim*m);
        }

        //  ''dim'' times of consecutive Householder Orthogonalizations
        for(j = i*dim; j < i*dim+dim; j++){
            //  vector h for Householder transformation
            memcpy(h, R + m*j + j, (m-j)*sizeof(double));
            int sign = (h[0] >= 0) - (h[0] < 0); 
            double v_nrm1, v_nrm2;
            la = m-j; ica = 1;  
            v_nrm1 = dnrm2(&la, h, &ica); 
            h[0] = h[0] - (-sign * v_nrm1);
            v_nrm2 = dnrm2(&la, h, &ica);

            //  normalization of h
            la = m - j; scala = 1.0/v_nrm2; ica = 1;
            dscal(&la, &scala, h, &ica);    
     
            //  elimination of the jth sub-column
            *(R + m*j + j) = - sign * v_nrm1;
            memset(R + m*j + j+1, 0, (m-j-1) * sizeof(double));   

            //  orthogonalization of columns right to the jth column
            trans = 'T';
            la = m-j; lb = n-j-1;
            scala = 1.0; scalb = 0.0; 
            lda = m; icb = 1; icc = 1;
            dgemv(&trans, &la, &lb, &scala, R + m*(j+1) + j, &lda, h, &icb, &scalb, mv, &icc);
            scala = -2.0; 
            dger(&la, &lb, &scala, h, &icb, mv, &icc, R + m*(j+1) + j, &lda);
        }

        //  Update the blcknorm
        if (dim == 3)
        {
            for(k = i + 1; k < n_blk; k++)
            {
                //  skip small column blocks
                if (blknorm[k] < coltol)
                {
                    blknorm[k] = 0;
                    continue; 
                }
                //  update column block norm
                double *init = R + m*dim*k + i*dim;  
                double tmp = blknorm[k]*blknorm[k] - 
                    (*(init + 0)) * (*(init + 0)) -
                    (*(init + 1)) * (*(init + 1)) - 
                    (*(init + 2)) * (*(init + 2)) - 
                    (*(init + m + 0)) * (*(init + m + 0)) -
                    (*(init + m + 1)) * (*(init + m + 1)) - 
                    (*(init + m + 2)) * (*(init + m + 2)) - 
                    (*(init + 2*m + 0)) * (*(init + 2*m + 0)) -
                    (*(init + 2*m + 1)) * (*(init + 2*m + 1)) - 
                    (*(init + 2*m + 2)) * (*(init + 2*m + 2));
                //  recalculation if blknorm[k] is too small
                if (tmp <= 1e-10)
                {
                    la = m-i*dim-dim; ica = 1;
                    double tmp1 = dnrm2(&la, R+m*(k*dim  )+i*dim+dim, &ica);
                    double tmp2 = dnrm2(&la, R+m*(k*dim+1)+i*dim+dim, &ica);
                    double tmp3 = dnrm2(&la, R+m*(k*dim+2)+i*dim+dim, &ica);
                    blknorm[k] = sqrt(tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3);
                }
                else
                    blknorm[k] = sqrt(tmp);                    
            }

        } else {
            for(k = i + 1; k < n_blk; k++)
            {
                /*  Need to modify  */
                //  skip small column blocks
                if (blknorm[k] < coltol)
                {
                    blknorm[k] = 0;
                    continue; 
                }
                //  update column block norms
                double tmp = 0.0;
                double tmp1 = 0.0;
                la = dim; ica = 1;
                for(l = k*dim; l < k*dim+dim; l++){
                    tmp1 = dnrm2(&la, R + m*l + i*dim, &ica);
                    tmp += tmp1 * tmp1;                    
                }
                blknorm[k] = sqrt(blknorm[k]*blknorm[k] - tmp);
                
            }
        }
        
        //  check full rank
        if (i == max_iter-1)
        {
            rank = max_iter*dim;
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

    /*  free space  */
    mxFree(blknorm);
    mxFree(h);
    mxFree(mv);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* check for proper number of arguments */
    if (nrhs != 4) {
        mexErrMsgTxt("QRdim requires 4 arguments.");
    }
    if (!mxIsDouble(prhs[0]) || mxIsSparse(prhs[0])) {
        mexErrMsgTxt( "Input must be a full matrix." );
    }
    if (!mxIsScalar(prhs[1])) {
        mexErrMsgTxt( "Input dimension must be an integer scalar.");
    }
    if (!mxIsChar(prhs[2])){
        mexErrMsgTxt( "Input stopping criteria must be a string as \"tol\" or \"rank\"." );
    }
    if (!mxIsScalar(prhs[3])){
        mexErrMsgTxt( "Input parameter for stopping must be a scalar.");
    }

    if (nlhs == 2 || nlhs == 3) {
        QRdim_without_Q(nlhs, plhs, nrhs, prhs);
    }
    // else if (nlhs == 3) {
    //     QRdim_with_Q(nlhs, plhs, nrhs, prhs);
    // }
    else {
        mexErrMsgTxt( "QR partial with Q component is not supported yet." );
    }
}
