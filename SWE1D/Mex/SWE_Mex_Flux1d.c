#include "mex.h"

#define real double

void SWE_NodalFlux1d(real hcrit, real gra,
                     real h, real qx,
                     real *Fh, real *Fq);

/* Computation of the flux term for two dimensional SWE.
 * Usages:
 * 	  [Fh, Fqx] = SWE_Flux1d(hmin, gra, h, qx)
 */
void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

	/* check input & output */
	if (nrhs != 4)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 2)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	real hmin = mxGetScalar(prhs[0]);
	real gra  = mxGetScalar(prhs[1]);
	real *h   = mxGetPr(prhs[2]);
	real *qx  = mxGetPr(prhs[3]);
	/* get dimensions */
	size_t np, ne;
	np = mxGetM(prhs[2]); 
	ne = mxGetN(prhs[2]);
	/* allocate output array */
	plhs[0] = mxCreateDoubleMatrix((mwSize)np, (mwSize)ne, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)np, (mwSize)ne, mxREAL);

	real *Fh  = mxGetPr(plhs[0]); 
    real *Fqx = mxGetPr(plhs[1]); 

	int i,j,ind=0;
	for (i=0;i<ne;i++){
		for(j=0;j<np;j++){
			SWE_NodalFlux1d(hmin, gra, h[ind], qx[ind],
                    Fh+ind, Fqx+ind);
			ind++;
		}
	}
    
    return;
}




void SWE_NodalFlux1d(real hcrit, real gra,
                     real h,   real qx,
                     real *Fh, real *Fq){
    if(h>hcrit){
        *Fh  = qx;
        *Fq = (real)(qx*qx/h + 0.5*gra*h*h);
    }else{
        *Fh  = 0; *Fq = 0;
    }

    return;
}