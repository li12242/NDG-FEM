#include "swe.h"

/* Computation of the flux term for two dimensional SWE.
 * Usages:
 * 	  [Fh, Fqx] = swe_nodal_flux(hmin, gra, h, qx, EToR)
 */
void mexFunction( int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[] )
{

    /* check input & output */
	if (nrhs != 5) mexErrMsgTxt("The number of input arguments should be 5.");
	if (nlhs != 2) mexErrMsgTxt("The number of output arguments should be 2.");

	/* get inputs */
	double hmin = mxGetScalar(prhs[0]);
	double gra  = mxGetScalar(prhs[1]);
	double *h   = mxGetPr(prhs[2]);
	double *qx  = mxGetPr(prhs[3]);
    signed char *etype = (signed char *)mxGetData(prhs[4]);
	/* get dimensions */
	size_t Np, K;
	Np = mxGetM(prhs[2]);
	K = mxGetN(prhs[2]);
	/* allocate output array */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);

	double *Fh  = mxGetPr(plhs[0]);
    double *Fq = mxGetPr(plhs[1]);

	int n,k,ind=0;
	for (k=0;k<K;k++){
        if (etype[k] == DRY){ // cell is dry
            for(n=0;n<Np;n++){
                Fh[ind] = 0.0;
                Fq[ind] = 0.0;
                ind ++;
            }
        }else{ // cell is wet
            for(n=0;n<Np;n++){
    			nodal_flux(hmin, gra, h[ind], qx[ind], Fh+ind, Fq+ind);
    			ind++;
    		}
        }
	}
    return;
}
