#include "mex.h"

#define real double

void SWE_NodalFlux2d(real hcrit, real gra,
                     real h, real qx, real qy,
                     real *Eh, real *Eqx, real *Eqy,
                     real *Gh, real *Gqx, real *Gqy);

/* Computation of the flux term for two dimensional SWE.
 * Usages:
 * 	  [Eh, Eqx, Eqy, Gh, Gqx, Gqy] = SWE_Flux2d(hmin, gra, h, qx, qy)
 */
void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

	// #define hmin prhs[0]
	// #define gra  prhs[1]
	// #define h    prhs[2]
	// #define qx   prhs[3]
	// #define qy   prhs[4]

	// #define Fh   plhs[0]
	// #define Fqx  plhs[1]
	// #define Fqy  plhs[2]
	// #define Gh   plhs[3]
	// #define Gqx  plhs[4]
	// #define Gqx  plhs[5]

	/* check input & output */
	if (nrhs != 5)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 6)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	real hmin = mxGetScalar(prhs[0]);
	real gra  = mxGetScalar(prhs[1]);
	real *h   = mxGetPr(prhs[2]);
	real *qx  = mxGetPr(prhs[3]);
	real *qy  = mxGetPr(prhs[4]);
	/* get dimensions */
	size_t np, ne;
	np = mxGetM(prhs[2]); 
	ne = mxGetN(prhs[2]);
	/* allocate output array */
	plhs[0] = mxCreateDoubleMatrix((mwSize)np, (mwSize)ne, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)np, (mwSize)ne, mxREAL);
	plhs[2] = mxCreateDoubleMatrix((mwSize)np, (mwSize)ne, mxREAL);
	plhs[3] = mxCreateDoubleMatrix((mwSize)np, (mwSize)ne, mxREAL);
	plhs[4] = mxCreateDoubleMatrix((mwSize)np, (mwSize)ne, mxREAL);
	plhs[5] = mxCreateDoubleMatrix((mwSize)np, (mwSize)ne, mxREAL);

	real *Eh  = mxGetPr(plhs[0]); 
    real *Eqx = mxGetPr(plhs[1]); 
    real *Eqy = mxGetPr(plhs[2]);
	real *Gh  = mxGetPr(plhs[3]); 
    real *Gqx = mxGetPr(plhs[4]); 
    real *Gqy = mxGetPr(plhs[5]);

	int i,j,ind=0;
	for (i=0;i<ne;i++){
		for(j=0;j<np;j++){
			SWE_NodalFlux2d(hmin, gra, h[ind], qx[ind], qy[ind], 
                    Eh+ind, Eqx+ind, Eqy+ind, Gh+ind, Gqx+ind, Gqy+ind);
			ind++;
		}
	}
    
    return;
}




void SWE_NodalFlux2d(real hcrit, real gra,
                     real h, real qx, real qy,
                     real *Eh, real *Eqx, real *Eqy,
                     real *Gh, real *Gqx, real *Gqy){

    if(h>hcrit){
        *Eh  = qx;
        *Eqx = (real)(qx*qx/h + 0.5*gra*h*h);
        *Eqy = qx*qy/h;
        *Gh  = qy;
        *Gqx = qx*qy/h;
        *Gqy = (real)(qy*qy/h + 0.5*gra*h*h);
    }else{
        *Eh  = 0; *Eqx = 0; *Eqy = 0;
        *Gh  = 0; *Gqx = 0; *Gqy = 0;
    }
}