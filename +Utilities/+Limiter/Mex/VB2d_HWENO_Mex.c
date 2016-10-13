#include "Limiter.h"

#define INFITY 10e10
#define TOTALERR 1e-12

/**
 * @brief
 * Use minmod function to limit the gradient and get the linear
 * limited result.
 *
 * Usages:
 *		shape = mesh.Shape;
 * 		hlim  = VB2d_Mex(h, mesh.J, shape.M, mesh.EToV, ...
 *			mesh.Shape.Fmask, mesh.Nv, mesh.x, mesh.y)
 *				
 */

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]){

	/* check input & output */
	if (nrhs != 8)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	real *h    = mxGetPr(prhs[0]);
	real *J    = mxGetPr(prhs[1]);
	real *M    = mxGetPr(prhs[2]);
	real *EToV = mxGetPr(prhs[3]);
	real *Fmask= mxGetPr(prhs[4]);
	real  Nv   = mxGetScalar(prhs[5]);
	real *x    = mxGetPr(prhs[6]);
	real *y    = mxGetPr(prhs[7]);

	int Nvert = (int) Nv;
	/* get dimensions */
	size_t Np, K;
	Np = mxGetM(prhs[0]);
	K  = mxGetN(prhs[0]);
	size_t Nfaces,Nfp;
	Nfaces = mxGetM(prhs[4]);
	Nfp    = mxGetN(prhs[4]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	real *hlim = mxGetPr(plhs[0]);

	/* cell averages */
	real *hmean = (real*) malloc(sizeof(real)*K );
	real *area  = (real*) malloc(sizeof(real)*K );
	real *xmean = (real*) malloc(sizeof(real)*K );
	real *ymean = (real*) malloc(sizeof(real)*K );
	//elemental integral coefficient
	real *w     = (real*) malloc(sizeof(real)*Np);
	/* volume integral coefficient */
	int i,j;
	for(i=0;i<Np;i++){
		w[i] = 0.0;
		for(j=0;j<Np;j++){
			w[i] += M[i*Np + j];
		}
	}

	// calculate volume mean value
	int k;
	for(k=0;k<K;k++){
		real *fld[3], cmean[3];
		fld[0] = h+k*Np; fld[1] = x+k*Np; fld[2] = y+k*Np;
		cellMean(Np, 3, fld, w, J+k*Np, cmean, area+k);
		hmean[k] = cmean[0];
		xmean[k] = cmean[1];
		ymean[k] = cmean[2];
	}

	// mexPrintf("Nv=%d\n", Nvert);
	/* vertex max/min bounds */
	real *hvmax = (real*) malloc(sizeof(real)*Nvert );
	real *hvmin = (real*) malloc(sizeof(real)*Nvert );
	// initialization
	int n;
	for(n=0;n<Nvert;n++){
		hvmax[n] = - INFITY;
		hvmin[n] = INFITY;
	}

	/* max and min value of each vertex */
	int f,v;
	for(k=0;k<K;k++){
		for(f=0;f<Nfaces;f++){
			i = (int) EToV[k + f*K]-1; // vertex index
			hvmax[i] = max( hvmax[i], hmean[k]);
			hvmin[i] = min( hvmin[i], hmean[k]);
			// mexPrintf("k=%d, v=%d, hmean=%f, max=%f, min=%f\n",
				// k, f, hmean[k],hvmax[i], hvmin[i]);
		}
	}

	// for(n=0;n<Nvert;n++){
	// 	mexPrintf("n=%d, hvmin=%f, hvmax=%f\n", n, hvmin[n], hvmax[n]);
	// }

	real alpha, dhdx, dhdy;
	real *hv = (real*) malloc(sizeof(real)*Nfaces );
	real *xv = (real*) malloc(sizeof(real)*Nfaces );
	real *yv = (real*) malloc(sizeof(real)*Nfaces );
	for(k=0;k<K;k++){
		alpha = 1.0;
		real hc = hmean[k];
		real xc = xmean[k];
		real yc = ymean[k];
		for(f=0;f<Nfaces;f++){
			i = k*Np + (int) Fmask[f]-1; // node index
			v = (int) EToV[k + f*K]-1; // vertex index
			xv[f] = x[i];
			yv[f] = y[i];
			hv[f] = h[i];
			if(hv[f]>hvmax[v]){
				hv[f]=hvmax[v];
			}else if(hv[f]<hvmin[v]){
				hv[f]=hvmin[v];
			}
		}

		HWENO_meanGradient(Nfaces, xv, yv, hv, xc, yc, hc, &dhdx, &dhdy);
		// mexPrintf("k=%d, dhdx=%f, dhdy=%f\n", k, dhdx, dhdy);
		getLocalVar(Np, hc, xc, yc, x+k*Np, y+k*Np, dhdx, dhdy, hlim+k*Np);

		// for(f=0;f<Nfaces;f++){
		// 	i = k*Np + (int) Fmask[f]-1; // node index
		// 	v = (int) EToV[k + f*K]-1; // vertex index
		// 	if (hlim[i] > hc*(1+TOTALERR) ){
		// 		alpha = min(alpha, ( hvmax[v]-hc )/(hlim[i] - hc) );
		// 	}else if(hlim[i] < hc*(1-TOTALERR)){
		// 		alpha = min(alpha, ( hvmin[v]-hc )/(hlim[i] - hc) );
		// 	}
		// 	// mexPrintf("k=%d, v=%d, ind=%d, hi=%f, hmax=%f, hmin=%f, alpha=%f\n",
		// 	// 	k, f, v, h[i], hvmax[v], hvmin[v], alpha);
		// }
		// // mexPrintf("k=%d, alpha=%f\n",k, alpha);
		// // get limited gradient
		// dhdx *= alpha;
		// dhdy *= alpha;
		// mexPrintf("k=%d, dhdx=%f, dhdy=%f\n", k, dhdx, dhdy);
		// getLocalVar(Np, hc, xc, yc, x+k*Np, y+k*Np, dhdx, dhdy, hlim+k*Np);
	}

	free(w);
	free(hv); free(xv); free(yv);

	free(hmean); free(area);
	free(xmean); free(ymean);
	free(hvmax); free(hvmin);
}
