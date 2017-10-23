#include "mex.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * brief: Calculate the mean value of each elements.
 * usages:
 *  [ c_mean ] = cell_mean( f_Q, w, J )
 */
void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    /* check input & output */
	if (nrhs != 3) mexErrMsgTxt("The number of input shoule be 3.");
	if (nlhs != 1) mexErrMsgTxt("The number of output shoule be 1.");

    /* get inputs */
	double *u = mxGetPr(prhs[0]);
	double *w = mxGetPr(prhs[1]);
	double *J = mxGetPr(prhs[2]);

	/* get dimensions */
    size_t Np = mxGetM(prhs[0]);
    size_t K  = mxGetN(prhs[0]);

	/* allocate output array */
	plhs[0] = mxCreateDoubleMatrix((mwSize)1, (mwSize)K, mxREAL);

	double *u_m = mxGetPr(plhs[0]);

    int n,k;
    #ifdef _OPENMP
    #pragma omp parallel for private(n) num_threads(DG_THREADS)
    #endif
    for(k=0;k<K;k++){
        u_m[k] = 0.0;
        double area = 0.0;
        for(n=0;n<Np;n++){
            int sk = k*Np + n;
            area += w[n]*J[sk];
            u_m[k] += w[n]*J[sk]*u[sk];
        }
        u_m[k] /= area;
    }
}
