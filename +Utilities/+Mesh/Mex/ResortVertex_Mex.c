#include "mex.h"
#include "VertexSort.h"

/*
 * Usages:
 * 		newEToV = ResortVertex_Mex(EToV, x, y);
 * 
 */


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]){


	/* check input & output */
	if (nrhs != 3)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	real *EToV  = mxGetPr(prhs[0]);
	real *x     = mxGetPr(prhs[1]);
	real *y     = mxGetPr(prhs[2]);

	/* get dimensions */
	size_t Nv, K;
	K  = mxGetM(prhs[0]);
	Nv = mxGetN(prhs[0]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)K, (mwSize)Nv, mxREAL);
	real *newEToV = mxGetPr(plhs[0]);

	int k;
	for(k=0;k<K;k++){
		VertexSort(K, Nv, EToV+k, x, y, newEToV+k);
	}

	return;
}