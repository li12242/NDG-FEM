#include "mex.h"
#include "polylib.h"

/*
* Returns value and derivative of Jacobi poly. at point z
* 
* This function returns the vectors poly_in and poly_d
*     containing the value of the $ n^th $ order Jacobi polynomial
*     $ P^{\alpha,\beta}_n(z) \alpha > -1, \beta > -1 $ and its
*     derivative at the np points in z[i]
* Input:
*       int         n   nth order Jacobi polynomial
*       double[np]  r   np points $\in$ [-1,1]
* Output:
* 		double[np] 	p
* 		double[np] 	dp
* Usages: [p,dp] = jacobfd(r, alpha, beta, n)
*       
*/ 

void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

	#define r_in prhs[0]
	#define alpha_in prhs[1]
	#define beta_in prhs[2]
	#define n_in prhs[3]

	#define p_out plhs[0]
	#define dp_out plhs[1]

	double *p, *dp, *r, alpha, beta;
	int n, np;

	/* check input & output */
	if (nrhs != 4)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 2)
		mexErrMsgTxt("Too many output arguments");

	np = (int)mxGetNumberOfElements(r_in);

	n = (int) mxGetScalar(n_in);
	alpha = mxGetScalar(alpha_in);
	beta = mxGetScalar(beta_in);
	r = mxGetPr(r_in);

	/* allocate output variable */
	p_out = mxCreateDoubleMatrix(np, 1, mxREAL);
	dp_out = mxCreateDoubleMatrix(np, 1, mxREAL);
	/* transfer Matlab variable to double */
	p = mxGetPr(p_out);
	dp = mxGetPr(dp_out);

	/* void jacobfd(int np, double *z, double *poly_in, double *polyd, int n, 
	     double alpha, double beta){ */

	jacobfd(np, r, p, dp, n, alpha, beta);

	return;
}