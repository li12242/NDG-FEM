#include "mex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define INFITY 10e10
#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

#define DEBUG 0
/*
 * Find the max and min vertex values through all the elements containg this vertex.
 * 
 * Usages:
 * 	[fmax, fmin] = vertex_extreme(Kv, VToE, c_mean)
 */
void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	/* check input & output */
	if (nrhs != 3)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 2)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	double *Kv = mxGetPr(prhs[0]);
	double *VToE = mxGetPr(prhs[1]);
	double *c_mean = mxGetPr(prhs[2]);

	size_t maxNe = mxGetM(prhs[1]);
	size_t Nv = mxGetN(prhs[1]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Nv, (mwSize)1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)Nv, (mwSize)1, mxREAL);
	double *f_max = mxGetPr(plhs[0]);
	double *f_min = mxGetPr(plhs[1]);

	int n, m;
	for (n = 0; n < Nv; n++)
	{ // initialization
		f_max[n] = -INFITY;
		f_min[n] = INFITY;
	}

#ifdef _OPENMP
#pragma omp parallel for private(m) num_threads(DG_THREADS)
#endif
	for (n = 0; n < Nv; n++)
	{
		int Ne = Kv[n]; // # of elemnts containing the vertex.
		double *vtoe = VToE + n * maxNe;

#if DEBUG
		mexPrintf("n=%d, Ne=%d\n", n, Ne);
#endif
		for (m = 0; m < Ne; m++)
		{
			int eid = (int)vtoe[m] - 1; // element index
			f_max[n] = max(f_max[n], c_mean[eid]);
			f_min[n] = min(f_min[n], c_mean[eid]);

#if DEBUG
			mexPrintf("m=%d, eid=%d, c_mean=%f\n", m, eid, c_mean[eid]);
#endif
		}
	}

	return;
}