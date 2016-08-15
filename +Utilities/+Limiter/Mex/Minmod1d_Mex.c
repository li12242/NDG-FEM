#include "Limiter.h"

/**
 * @brief
 * Use minmod function to limit the gradient and get the linear 
 * limited result.
 * 
 * Usages:
 *		shape = mesh.Shape;
 * 		hlim  = Minmod1d_Mex...
 *			(h, mesh.J, shape.M, shape.Fmask, mesh.EToE, mesh.x)
 */
void mexFunction (int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

	/* check input & output */
	if (nrhs != 6)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	real *h    = mxGetPr(prhs[0]);
	real *J    = mxGetPr(prhs[1]);
	real *M    = mxGetPr(prhs[2]);
	real *Fmask= mxGetPr(prhs[3]);
	real *EToE = mxGetPr(prhs[4]);
	real *x    = mxGetPr(prhs[5]);
	
	/* get dimensions */
	size_t Np, K;
	Np = mxGetM(prhs[0]); 
	K  = mxGetN(prhs[0]);
	size_t Nfaces,Nfp;
	Nfaces = mxGetM(prhs[3]);
	Nfp    = mxGetN(prhs[3]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	real *hlim = mxGetPr(plhs[0]);

	/* cell averages */
	real *hmean = (real*) malloc(sizeof(real)*K );
	real *area  = (real*) malloc(sizeof(real)*K );
	//elemental integral coefficient
	real *w     = (real*) malloc(sizeof(real)*Np); 
	int i,j,k,sk;
	for(i=0;i<Np;i++){
		w[i] = 0.0;
		for(j=0;j<Np;j++){
			w[i] += M[i*Np + j];
		}
	}
	for (k=0;k<K;k++){
		area [k] = 0.0;
		hmean[k] = 0.0;
		for(i=0;i<Np;i++){
			sk = (k*Np + i);
			hmean[k] += w[i]*J[sk]*h[sk];
			area [k] += w[i]*J[sk];
		}
		hmean[k] /= area[k];
	}

	// reconstruct over each cell
	int  v1  = (int) Fmask[0] - 1;
	int  v2  = (int) Fmask[1] - 1;
	for(k=0;k<K;k++){
		/* cell gradient */
		int sk1 = k*Np + v1;
		int sk2 = k*Np + v2;
		real xc  =(x[sk2] + x[sk1])/2.0;
		real dh  = h[sk2] - h[sk1];
		/* limited gradient */
		int e1 = (int)EToE[k  ] - 1;
		int e2 = (int)EToE[k+K] - 1;
		real a[3];
		a[0] = dh/area[k];
		a[1] = (hmean[e2] - hmean[k] )/area[k];
		a[2] = (hmean[k]  - hmean[e1])/area[k];
		real dhlim;
		minmod(3, a, &dhlim);
		/* reconstruction */
		for(i=0;i<Np;i++){
			sk = (k*Np + i);
            real delta = (real)(x[sk] - xc);
            hlim[sk] = hmean[k] + delta*dhlim;
		}
	}

	free(w); free(hmean); free(area);
	return;
}