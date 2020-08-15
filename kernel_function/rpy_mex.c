#include <math.h>
#include "mex.h"
#include "matrix.h"

#define M_PI 3.14159265358979323846

void rpy(int np, const double *pos, double a, double eta, double *A)
{
    int i;
    const double C = 1./(6 * M_PI * a * eta) ; 
    for (i=0; i < np; i++)
    {
        double posi[4];
        posi[0] = pos[i  ];
        posi[1] = pos[i + np];
        posi[2] = pos[i + 2*np];
        
        const int ld = 3*np;
        int base; 
#define mat(k, l) A[base + ld*k + l]
        base = 3*i*ld + 3*i;
        mat(0,0) = C;
        mat(0,1) = 0;
        mat(0,2) = 0;
        mat(1,0) = 0;
        mat(1,1) = C;
        mat(1,2) = 0;
        mat(2,0) = 0;
        mat(2,1) = 0;
        mat(2,2) = C;
        
        int j;
        for (j=i+1; j<np; j++)
        {
            double rvec[4];
            double s, s2;
            
            rvec[0] = posi[0] - pos[j     ];
            rvec[1] = posi[1] - pos[j+np  ];
            rvec[2] = posi[2] - pos[j+2*np];            
            s2 = rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2];
            s = sqrt(s2);
            rvec[0] /= s;
            rvec[1] /= s;
            rvec[2] /= s;
            
            double t1, t2;
            if (s < 2 * a)
            {
                t1 = C * (1- 9./32.*s/a);
                t2 = C * 3. / 32. * s / a;    
            }
            else 
            {
                t1 = C * 3./4.*a/s * (1 + 2./3.*a*a/s2);
                t2 = C * 3./4.*a/s * (1 - 2.*a*a/s2); 
            }
            
            base = 3*i*ld + 3*j;
            mat(0,0) = t2*rvec[0]*rvec[0] + t1;
            mat(0,1) = t2*rvec[0]*rvec[1];
            mat(0,2) = t2*rvec[0]*rvec[2];
            mat(1,0) = t2*rvec[1]*rvec[0];
            mat(1,1) = t2*rvec[1]*rvec[1] + t1;
            mat(1,2) = t2*rvec[1]*rvec[2];
            mat(2,0) = t2*rvec[2]*rvec[0];
            mat(2,1) = t2*rvec[2]*rvec[1];
            mat(2,2) = t2*rvec[2]*rvec[2] + t1;

            base = 3*j*ld + 3*i;
            mat(0,0) = t2*rvec[0]*rvec[0] + t1;
            mat(0,1) = t2*rvec[0]*rvec[1];
            mat(0,2) = t2*rvec[0]*rvec[2];
            mat(1,0) = t2*rvec[1]*rvec[0];
            mat(1,1) = t2*rvec[1]*rvec[1] + t1;
            mat(1,2) = t2*rvec[1]*rvec[2];
            mat(2,0) = t2*rvec[2]*rvec[0];
            mat(2,1) = t2*rvec[2]*rvec[1];
            mat(2,2) = t2*rvec[2]*rvec[2] + t1;
        }
    }
    return ;
}

void rpy_off_block(int np1, const double *pos1, int np2, const double *pos2, double a, double eta, double *A)
{
    int i;
    const double C = 1./(6 * M_PI * a * eta) ; 
    for (i=0; i < np1; i++)
    {
        double posi[4];
        posi[0] = pos1[i  ];
        posi[1] = pos1[i+np1];
        posi[2] = pos1[i+2*np1];
        
        const int ld = 3*np1;
        int base; 
#define mat(l, k) A[base + ld*k + l] 
        
        int j;
        for (j=0; j<np2; j++)
        {
            double rvec[4];
            double s, s2;            
          
            rvec[0] = posi[0] - pos2[j  ];
            rvec[1] = posi[1] - pos2[j+np2];
            rvec[2] = posi[2] - pos2[j+2*np2];           
            s2 = rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2];
            s = sqrt(s2);
            base = 3*j*ld + 3*i;
            
            if (s < 1e-10)
            {
                mat(0,0) = C;
                mat(0,1) = 0;
                mat(0,2) = 0;
                mat(1,0) = 0;
                mat(1,1) = C;
                mat(1,2) = 0;
                mat(2,0) = 0;
                mat(2,1) = 0;
                mat(2,2) = C;            
                continue;
            }
                
            rvec[0] /= s;
            rvec[1] /= s;
            rvec[2] /= s;
 
            double t1, t2;
            if (s < 2 * a)
            {
                t1 = C * (1- 9./32.*s/a);
                t2 = C * 3. / 32. * s / a;    
            }
            else 
            {
                t1 = C * 3./4./s * (1 + 2./3.*a*a/s2);
                t2 = C * 3./4./s * (1 - 2.*a*a/s2); 
            }
            
            mat(0,0) = t2*rvec[0]*rvec[0] + t1;
            mat(0,1) = t2*rvec[0]*rvec[1];
            mat(0,2) = t2*rvec[0]*rvec[2];
            mat(1,0) = t2*rvec[1]*rvec[0];
            mat(1,1) = t2*rvec[1]*rvec[1] + t1;
            mat(1,2) = t2*rvec[1]*rvec[2];
            mat(2,0) = t2*rvec[2]*rvec[0];
            mat(2,1) = t2*rvec[2]*rvec[1];
            mat(2,2) = t2*rvec[2]*rvec[2] + t1;
        }
    }
    return ;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int np1,np2;
    const double *pos1;
    const double *pos2;
    double a, eta;
    double *A;
    
    if (mxIsCell(prhs[0]))
    {
        np1 = mxGetM(mxGetCell(prhs[0], 0));
        pos1 = mxGetPr(mxGetCell(prhs[0], 0));
        np2 = mxGetM(mxGetCell(prhs[0], 1));
        pos2 = mxGetPr(mxGetCell(prhs[0], 1));
        a = mxGetScalar(prhs[1]);
        eta = mxGetScalar(prhs[2]);
        
        // allocate space for output
        plhs[0] = mxCreateDoubleMatrix(np1*3, np2*3, mxREAL);
        A = mxGetPr(plhs[0]);
        
        rpy_off_block(np1, pos1, np2, pos2, a, eta, A);
    }
    else
    {
        //read the input
        np1 = mxGetM(prhs[0]);
        pos1 = mxGetPr(prhs[0]);
        a = mxGetScalar(prhs[1]);
        eta = mxGetScalar(prhs[2]);
        
        // allocate space for output
        plhs[0] = mxCreateDoubleMatrix(np1*3, np1*3, mxREAL);
        A = mxGetPr(plhs[0]);
        rpy(np1, pos1, a, eta, A);
    }
}
