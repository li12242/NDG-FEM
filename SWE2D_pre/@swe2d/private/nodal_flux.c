#include "swe.h"

/* Computation of the flux term for two dimensional SWE.
 * Usages:
 * 	  [Eh, Eqx, Eqy, Gh, Gqx, Gqy] = swe_nodal_flux(hmin, gra, h, qx, qy, z, EToR)
 */
void mexFunction( int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[] )
{

    /* check input & output */
	if (nrhs != 7) mexErrMsgTxt("The number of input arguments should be 5.");
	if (nlhs != 6) mexErrMsgTxt("The number of output arguments should be 2.");

	/* get inputs */
	double hmin = mxGetScalar(prhs[0]);
	double gra  = mxGetScalar(prhs[1]);
	double *h   = mxGetPr(prhs[2]);
	double *qx  = mxGetPr(prhs[3]);
    double *qy = mxGetPr(prhs[4]);
    double *z = mxGetPr(prhs[5]);
    signed char *etype = (signed char *)mxGetData(prhs[6]);
	/* get dimensions */
	size_t Np, K;
	Np = mxGetM(prhs[2]);
	K = mxGetN(prhs[2]);
	/* allocate output array */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
    plhs[2] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	plhs[3] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
    plhs[4] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	plhs[5] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);

	double *Eh  = mxGetPr(plhs[0]);
    double *Eqx = mxGetPr(plhs[1]);
    double *Eqy = mxGetPr(plhs[2]);
    double *Gh  = mxGetPr(plhs[3]);
    double *Gqx = mxGetPr(plhs[4]);
    double *Gqy = mxGetPr(plhs[5]);

	int n,k;
    #ifdef _OPENMP
    #pragma omp parallel for private(n) num_threads(DG_THREADS)
    #endif
	for (k=0;k<K;k++){
        int ind = k*Np;
        if (etype[k] == DRY){ // cell is dry
            for(n=0;n<Np;n++){
                Eh[ind] = 0.0; Eqx[ind] = 0.0; Eqy[ind] = 0.0;
                Gh[ind] = 0.0; Gqx[ind] = 0.0; Gqy[ind] = 0.0;
                ind ++;
            }
        }else{ // cell is wet
            for(n=0;n<Np;n++){
    			nodal_flux(hmin, gra, h[ind], qx[ind], qy[ind], z[ind],
                    Eh+ind, Eqx+ind, Eqy+ind,
                    Gh+ind, Gqx+ind, Gqy+ind);
    			ind++;
    		}
        }
	}
    return;
}
