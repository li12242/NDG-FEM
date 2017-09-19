#include "mex.h"
#include "vert_sort.h"

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
	Nv = mxGetM(prhs[0]);
	K  = mxGetN(prhs[0]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Nv, (mwSize)K, mxREAL);
	real *newEToV = mxGetPr(plhs[0]);

	int k;
	for(k=0;k<K;k++){
		VertexSort(Nv, EToV+k*Nv, x, y, newEToV+k*Nv);
	}

	return;
}