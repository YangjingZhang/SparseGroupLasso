/***********************************************************************
* Compute fz = mexfz(z,ind,grpNUM)
* by Xudong Li, 8 May 2017
***********************************************************************/

#include <math.h>
#include <mex.h>
#include <matrix.h>
#include <string.h> /* needed for memcpy() */


#if !defined(SQR)
#define SQR(x) ((x)*(x))
#endif


/**********************************************************
* 
***********************************************************/
double normfun(const double *z, int kstart, int kend)

{ 
    double nrm, nrm2 =0.0;
    ptrdiff_t i;
    
    for (i=kstart-1; i< kend; i++) {
        nrm2 += SQR(z[i]);
    }
    nrm = sqrt(nrm2);
    /*printf(" kstart = %d, kend = %d ", kstart, kend);*/
    return nrm;
}


void mexFunction(
      int nlhs,   mxArray  *plhs[], 
      int nrhs,   const mxArray  *prhs[] )

{        
    double   *ind, *z, *fz;
    double   nrm, tmp=0.0; 

    ptrdiff_t  grpNUM;
    int   m, j, k,kstart, kend; 

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

   if (nrhs > 3){
      mexErrMsgTxt("mexFnorm: requires at most 4 input arguments."); }
   if (nlhs > 1){ 
      mexErrMsgTxt("mexFnorm: requires at most 1 output argument."); }   

/* CHECK THE DIMENSIONS */
    
    
    m = mxGetM(prhs[0]); 
    z = mxGetPr(prhs[0]);
    ind = mxGetPr(prhs[1]);
    grpNUM = (int)*mxGetPr(prhs[2]);
    
    /***** create return argument *****/    
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); 
    fz = mxGetPr(plhs[0]); 
    

    /***** Do the computations in a subroutine *****/  
    for (j=0; j<grpNUM;j++) {
        kstart = (ptrdiff_t) ind[3*j];
        kend = (ptrdiff_t) ind[1+3*j];
        nrm = normfun(z, kstart, kend);
        tmp += ind[2+3*j]*nrm;     
    }
    fz[0] = tmp;
    return;
}
/**********************************************************/
