#include "VectorOperator.h"

/*
 * Usages:
 *   locFlag = FindLocCell_Mex(x, y, locVertList, xp, yp);
 *				
 */

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]){

	/* check input & output */
	if (nrhs != 5)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	double *x    = mxGetPr(prhs[0]);
	double *y    = mxGetPr(prhs[1]);
	double *locVertList = mxGetPr(prhs[2]);
	double *xp   = mxGetPr(prhs[3]);
	double *yp   = mxGetPr(prhs[4]);

	size_t Np, K;
	Np = mxGetM(prhs[0]);
	K  = mxGetN(prhs[0]);

	size_t Nvert;
	Nvert = mxGetM(prhs[2]);

	size_t Mp; // No. of unknown points
	Mp = mxGetM(prhs[3]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)1, (mwSize)K, mxREAL);
	double *isInside = mxGetPr(plhs[0]);

	if(Nvert<2)
		mexErrMsgTxt("Number of vertex is less than 2.");

	int i,m,k;
	for(m=0;m<Mp;m++){
		POINT p;
		p.x = xp[m]; p.y = yp[m];

		for(k=0;k<K;k++){
			double f = 1;
			for(i=0;i<Nvert;i++){
				int v  = (int) locVertList[i]-1;
				int v1 = (int) locVertList[(i-1+Nvert)%Nvert]-1;
				int v2 = (int) locVertList[(i+1)%Nvert]-1;

				POINT A,B,C;
				A.x = x[v +k*Np]; A.y = y[v +k*Np];
				B.x = x[v1+k*Np]; B.y = y[v1+k*Np];
				C.x = x[v2+k*Np]; C.y = y[v2+k*Np];

				POINT AB, AC, AP;
				minus(B, A, &AB);
				minus(C, A, &AC);
				minus(p, A, &AP);

				f *= cross(AB, AP)*cross(AB, AC);

				// mexPrintf("k=%d, Nv=%d, AB=[%f,%f], AC=[%f,%f], AP=[%f,%f], f=%f\n", 
				// 	k, i, AB.x, AB.y, AC.x, AC.y, AP.x, AP.y, f);

				if(f<0) break;
			}
			if(f>0) isInside[k] = m+1;
		}
	}

	return;
}