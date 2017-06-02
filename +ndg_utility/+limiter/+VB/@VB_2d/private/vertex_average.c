#include "mex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define INFITY 10e10
#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

#define DEBUG 0

void weight_vertex_value(int Nvert, int MaxKv, double *Kv,
    double *VToE, double *VToC, double *hc, double *hv)
{
    int n,k;
    #pragma omp parallel for private(k)
    for(n=0;n<Nvert;n++){ // loop for all vertex
        hv[n] = 0.0;
        int Ne = (int) Kv[n];

        double *vtoe = VToE + n*MaxKv; // nth column
        double *vtoc = VToC + n*MaxKv;
        for(k=0;k<Ne;k++){ // loop for adjacent elements
            int eind = (int) vtoe[k]-1; // to C type
            double w = vtoc[k];
            hv[n] += w*hc[eind];
        }
    }
    return;
}
/*
 * Find max and min vertex values.
 *
 * Usages:
 * 	[fv] = vertex_extreme(fc, Kv, VToE, VToC)
 */
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
    /* check input & output */
	if (nrhs != 4) mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1) mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	double *fc = mxGetPr(prhs[0]);
	double *Kv = mxGetPr(prhs[1]);
	double *VToE = mxGetPr(prhs[2]);
	double *VToC = mxGetPr(prhs[3]);

	/* get dimensions */
    size_t MaxKv = mxGetM(prhs[2]);
    size_t Nvert = mxGetN(prhs[2]);

    /* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Nvert, (mwSize)1, mxREAL);
	double *fv = mxGetPr(plhs[0]);
    weight_vertex_value(Nvert, MaxKv, Kv, VToE, VToC, fc, fv);

    return;
}
