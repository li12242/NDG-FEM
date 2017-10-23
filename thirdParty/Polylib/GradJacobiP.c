#include "mex.h"
#include "polylib.h"

/*
* Purpose: Evaluate the derivative of the Jacobi polynomial 
* 	of type (alpha,beta)>-1, at points r for order N and 
* 	returns dP[1:length(r))]
* 
* Input:
*       double[np]  r   np points $\in$ [-1,1]
* 		double 		alpha
* 		double 		beta
* 		int 		n
* Output:
* 		double[np] 	dp
*
* Usages: dp = GradJacobiP(r, alpha, beta, n)
*       
*/ 


void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

	#define r_in prhs[0]
	#define alpha_in prhs[1]
	#define beta_in prhs[2]
	#define n_in prhs[3]

	#define dp_out plhs[0]

	double *dp, *r, alpha, beta;
	int np, n;

	/* check input & output */
	if (nrhs != 4)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Too many output arguments");

	np = (int)mxGetNumberOfElements(r_in);

	n = (int) mxGetScalar(n_in);
	alpha = mxGetScalar(alpha_in);
	beta = mxGetScalar(beta_in);
	r = mxGetPr(r_in);

	/* allocate output variable */
	dp_out = mxCreateDoubleMatrix(np, 1, mxREAL);

	dp = mxGetPr(dp_out);

	/* void GradjacobiP(int np, double *z, double *dp, 
	int n, double alpha, double beta) */
	GradjacobiP(np, r, dp, n, alpha, beta);

	return;
}