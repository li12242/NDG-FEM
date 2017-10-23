#include "mex.h"
#include "polylib.h"

/* 
 * Get Gauss-Legendre points and weights.
 * Input: 	
 * 		np 	- number of nodes in [-1, 1]
 * Output:
 * 		z 	- x coordinate
 * 		w 	- weights
 * 		
 * Usages: [z,w] = zwglj(np)
*/

void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

	#define np_in prhs[0]
	#define z_out plhs[0]
	#define w_out plhs[1]

	double *z, *w;
	int np;

	/* check input & output */
	if (nrhs != 1) mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 2) mexErrMsgTxt("Too many output arguments");

	np = (int) mxGetScalar(np_in);
	z_out = mxCreateDoubleMatrix(np,1,mxREAL);
	w_out = mxCreateDoubleMatrix(np,1,mxREAL);

	z = mxGetPr(z_out);
	w = mxGetPr(w_out);

	/* void zwglj(double *z, double *w, int np, double alpha, double beta) */

	zwgl(z, w, np);

}