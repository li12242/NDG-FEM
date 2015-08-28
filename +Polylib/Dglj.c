#include "mex.h"
#include "polylib.h"

/* 
 * Compute the Derivative Matrix and its transpose associated with the 
 * 	Gauss-Lobatto-Jacobi zeros.
 *
 * Compute the derivative matrix, associated with the n_th order 
 * 	Lagrangian interpolants through the np Gauss-Lobatto-Jacobi 
 *  points z such that
 *  $\frac{du}{dz}(z[i]) =  \sum_{j=0}^{np-1} D[i*np+j] u(z[j])$
 *  $D[i*np+j] = \frac{\partial l_j}{\partial z} \right|_{z=z_i}$
 *
 * Input: 	
 * 		r 	- x coordinate
 * Output:
 * 		D 	- derivative matrix
 * 		
 * Usages: D = Dglj(r)
*/

void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

	#define r_in prhs[0]
	#define D_out plhs[0]

	double *D, *r, *Dt;
	int np;


	/* check input & output */
	if (nrhs != 1)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Too many output arguments");

	np = (int)mxGetNumberOfElements(r_in);
	r = mxGetPr(r_in);
	/* allocate output argument */
	D_out = mxCreateDoubleMatrix(np, np, mxREAL);
	D = mxGetPr(D_out);

	Dt = (double*)malloc(sizeof(double)*np*np);
	/* void Dglj(double *D, double *Dt, double *z, int np,
	  double alpha, double beta) */
	Dglj(D, Dt, r, np, 0, 0);

	free(Dt);
}