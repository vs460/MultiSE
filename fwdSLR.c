#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>  
#include "mex.h"
#include "math.h"
#include "matrix.h"
#define TWOPI	6.283185


void fwdSLR(double *rf_re,double *rf_im, int N, double *x, int Nf, double gamma, double dt, double *a_re,double *a_im, double *b_re,double *b_im)
{
/* Performs fwd SLR-transform, i.e. calculates the complex alpha and beta
 polynomials from the RF pulse

 INPUT:
   - rf_re + I*re_im = RF pulse
   - N               = length of the RF pulse
   - x               = position/frequency vector where a and b will be evaluated
   - Nf              = length of x
   - gamma           = gyromagnetic ratio
   - dt              = sampling time of the pulse
 OUTPUT:
   - a_re + I*a_im   = alpha polynomial
   - b_re + I*b_im   = beta polynomial 
 */

int k,j;
double _Complex alpha,beta,a_j,b_j;
double _Complex product[2];
double phi,C,nx,ny,nz,xx;

/* calculate a and b for each location */
for(k = 0; k < Nf; k++){
    xx = x[k];
    
    /* initialize alpha and beta */
    alpha = 1;
    beta = 0;
    
    /* loop to find a_j from a_j-1 */
    for(j = 0; j < N; j++){
        /* calculate auxiliary values for the algorithm */
        phi = -TWOPI*gamma*dt*sqrt(cabs(rf_re[j]+I*rf_im[j])*cabs(rf_re[j]+I*rf_im[j]) + (xx/gamma)*(xx/gamma));
        C  = TWOPI*gamma*dt/fabs(phi);
        nx = C*rf_re[j];
        ny = C*rf_im[j];    
        nz = C*xx/gamma;
        
        /* calculate rotation matrices */
        a_j = cos(phi/2) + I*nz*sin(phi/2);
        b_j = I*(nx+I*ny)*sin(phi/2);
        
        product[0] = a_j*alpha - conj(b_j)*beta;
        product[1] = b_j*alpha + conj(a_j)*beta;
        
        /* iterate for the next alpha and beta */
        alpha = product[0];
        beta  = product[1];
    }
    
    a_re[k] = creal(alpha);
    a_im[k] = cimag(alpha);
    b_re[k] = creal(beta); 
    b_im[k] = cimag(beta); 
}
}




void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    double *rf_re,*rf_im,*x;
    double *a_re,*a_im,*b_re,*b_im;
    double gamma,dt;
    int N,Nf;
      
    /*pointers to input arrays*/
    rf_re = mxGetPr(prhs[0]);
    rf_im = mxGetPr(prhs[1]);
    N     = mxGetScalar(prhs[2]);
    x     = mxGetPr(prhs[3]);
    Nf    = mxGetScalar(prhs[4]);
    gamma = mxGetScalar(prhs[5]);
    dt    = mxGetScalar(prhs[6]);
    
    /* auxiliary calculation and memory allocation for output */
    a_re = mxCalloc(Nf,sizeof(double));
    a_im = mxCalloc(Nf,sizeof(double));
    b_re = mxCalloc(Nf,sizeof(double));
    b_im = mxCalloc(Nf,sizeof(double));
         
    /* create pointer to output argument*/
    plhs[0] = mxCreateDoubleMatrix(1,Nf,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,Nf,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,Nf,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,Nf,mxREAL);
    
    /*pointers to output arrays*/
    a_re = mxGetPr(plhs[0]);
    a_im = mxGetPr(plhs[1]);
    b_re = mxGetPr(plhs[2]);
    b_im = mxGetPr(plhs[3]);
                 
    fwdSLR(rf_re,rf_im,N,x,Nf,gamma,dt,a_re,a_im,b_re,b_im);  
    
}
