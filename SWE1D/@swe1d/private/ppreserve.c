#include "swe.h"

/*
 * Usages:
 *		[h_pos, q_pos] = ppreserve(hcrit, h, q, hc, qc);
 */
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{

	/* check input & output */
	if (nrhs != 5)
		mexErrMsgTxt("The number of input shoule be 5.");
	if (nlhs != 2)
		mexErrMsgTxt("The number of output shoule be 2.");

	/* get inputs */
    double hcrit = mxGetScalar(prhs[0]);
	double *h = mxGetPr(prhs[1]);
	double *q = mxGetPr(prhs[2]);
	double *hc = mxGetPr(prhs[3]);
	double *qc = mxGetPr(prhs[4]);

	/* get dimensions */
	size_t Np = mxGetM(prhs[1]);
	size_t K  = mxGetN(prhs[1]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);

	double *h_pos = mxGetPr(plhs[0]);
	double *q_pos = mxGetPr(plhs[1]);
	double ksi = 0.0;
	int k,n,sk,ind;
	// cell area and scalar averages
	sk = 0;
	double hmean, qmean, theta;
	for (k = 0;k < K; k++){
		hmean = hc[k];
		qmean = qc[k];

		if(hmean<=ksi){
			for(n = 0; n < Np; n++ ){
				ind = k * Np + n;
				h[ind] = 0;
				q[ind] = 0;
			}
            continue;
		}

		double hmin = h[k*Np];
		for(n=0;n<Np;n++){ hmin = min(hmin, h[k*Np + n]); }

		if(hmin < hmean){
			theta = min( (hmean-ksi)/(hmean-hmin), 1.0 );
		}else{ theta = 0.0; qmean = 0.0; }

		for(n=0;n<Np;n++){
			ind = k*Np + n;
			h_pos[ind] = theta*(h[ind] - hmean) + hmean;
			q_pos[ind] = theta*(q[ind] - qmean) + qmean;

		    if(h_pos[ind] < hcrit){ q_pos[ind] = 0.0; } // dry nodes
		}

	}
	return;
}
