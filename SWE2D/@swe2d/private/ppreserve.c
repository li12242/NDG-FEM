#include "swe.h"

/*
 * Usages:
 *		[h_pos, qx_pos, qy_pos] = ppreserve(hcrit, h, qx, qy, hc, qxc, qyc);
 */
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{

	/* check input & output */
	if (nrhs != 7) mexErrMsgTxt("The number of input shoule be 5.");
	if (nlhs != 3) mexErrMsgTxt("The number of output shoule be 2.");

	/* get inputs */
    double hcrit = mxGetScalar(prhs[0]);
	double *h = mxGetPr(prhs[1]);
	double *qx = mxGetPr(prhs[2]);
	double *qy = mxGetPr(prhs[3]);
	double *hc = mxGetPr(prhs[4]);
	double *qxc = mxGetPr(prhs[5]);
	double *qyc = mxGetPr(prhs[6]);

	/* get dimensions */
	size_t Np = mxGetM(prhs[1]);
	size_t K  = mxGetN(prhs[1]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	plhs[2] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);

	double *h_pos = mxGetPr(plhs[0]);
	double *qx_pos = mxGetPr(plhs[1]);
	double *qy_pos = mxGetPr(plhs[2]);
	int k,n;
    double ksi = 0.0;
	// cell area and scalar averages
    #pragma omp parallel for private(n, ksi)
	for (k = 0;k < K; k++){
		double hmean = hc[k];
		double qxmean = qxc[k];
		double qymean = qyc[k];
		if(hmean<=ksi){
			for(n = 0; n < Np; n++ ){
				int ind = k * Np + n;
				h[ind] = 0;
				qx[ind] = 0;
				qy[ind] = 0;
			}
            continue;
		}

		double hmin = h[k*Np];
		for(n=0;n<Np;n++){ hmin = min(hmin, h[k*Np + n]); }

        double theta;
		if(hmin < hmean){
			theta = min( (hmean-ksi)/(hmean-hmin), 1.0 );
		}else{ theta = 0.0; }

		for(n=0;n<Np;n++){
			int ind = k*Np + n;
			h_pos[ind] = theta*(h[ind] - hmean) + hmean;
			qx_pos[ind] = theta*(qx[ind] - qxmean) + qxmean;
			qy_pos[ind] = theta*(qy[ind] - qymean) + qymean;

		    if(h_pos[ind] < hcrit){ // dry nodes
				qx_pos[ind] = 0.0; qy_pos[ind] = 0.0;
			}
		}

	}
	return;
}
