#include "mex.h"
#include "polylib.h"

/* 
 * get normalized Jacobi polynamial
 * Input: 	r 		- node coordinate in [-1, 1]
 * 		alpha
 * 		beta 	
 * 		n
 * Usages: p = JacobiP(r,alpha,beta,n)
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	#define r_in prhs[0]
	#define alpha_in prhs[1]
	#define beta_in prhs[2]
	#define n_in prhs[3]

	#define p_out plhs[0]

	int np, n;
	double alpha, beta;
	double *r, *poly;

	/* check input & output */
	if (nrhs != 4)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Too many output arguments");

	np = (int)mxGetNumberOfElements(r_in);
	n = (int)mxGetScalar(n_in);
	r = mxGetPr(r_in);
	alpha = mxGetScalar(alpha_in);
	beta = mxGetScalar(beta_in);

	/* allocate output argument */
	p_out = mxCreateDoubleMatrix(np, 1, mxREAL);
	poly = mxGetPr(p_out);

	/* jacobiP(int np, double *z, double *poly_in, int n,  d
	ouble alpha, double beta) */

	jacobiP(np, r, poly, n, alpha, beta);

}