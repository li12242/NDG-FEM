#include "mex.h"

#define real double

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

/*
 * 
 * Usages:
 *		shape  = mesh.Shape;
 *		[h, q] = SWE_Mex_PositivePreserving1d(h, q, shape.M, mesh.J, phys.minDepth);
 */
void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

	/* check input & output */
	if (nrhs != 5)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 2)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	real *h    = mxGetPr(prhs[0]);
	real *qx   = mxGetPr(prhs[1]);
	real *M    = mxGetPr(prhs[2]);
	real *J    = mxGetPr(prhs[3]);
	real hcri  = mxGetScalar(prhs[4]);

	/* get dimensions */
	size_t Np, K;
	Np = mxGetM(prhs[0]); 
	K  = mxGetN(prhs[0]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);

	real *hp  = mxGetPr(plhs[0]);
	real *qxp = mxGetPr(plhs[1]);

	real ksi = 0.0;

	//elemental integral coefficient
	real *w     = (real*) malloc(sizeof(real)*Np); 
	int i,j,k,sk,ind;
	for(i=0;i<Np;i++){
		w[i] = 0.0;
		for(j=0;j<Np;j++){
			w[i] += M[i*Np + j];
		}
	}

	// cell area and scalar averages
	sk = 0;
	real area, hmean, qxmean;
	real hmin, theta;
	for (i=0;i<K;i++){
		area  = 0.0;
		hmean = 0.0;
		qxmean = 0.0;
		for(j=0;j<Np;j++){
			hmean  += w[j]*J[sk]*h[sk];
			qxmean += w[j]*J[sk]*qx[sk];
 			area   += w[j]*J[sk];
			sk++;
		}
		hmean /= area;
		qxmean /= area;


		if(hmean<=ksi){
			for(j=0;j<Np;j++){
				ind = i*Np + j;
				h[ind] += (ksi - hmean);
				qx[ind] = 0;
			}
			hmean = ksi;
			qxmean = 0.0;
		}

		hmin = h[i*Np];
		for(j=0;j<Np;j++){
			ind = i*Np + j;
			hmin = min(hmin, h[ind]);
		}

		if(hmin < hmean){
			theta = min( (hmean-ksi)/(hmean-hmin), 1.0 );
		}else{
			theta = 0.0;
		}

		for(j=0;j<Np;j++){
			ind = i*Np + j;
			hp [ind] = theta*(h [ind] - hmean) + hmean;
			qxp[ind] = theta*(qx[ind] - qxmean) + qxmean;

		    if(hp[ind] < hcri){ // dry nodes
		    	qxp[ind] = 0.0;
		    }
		}

	}

	free(w);

	return;
}