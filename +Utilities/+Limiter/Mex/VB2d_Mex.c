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
 * 		hlim  = VB2d_Mex...
 *			(h, mesh.J, shape.M, mesh.EToV, mesh.Shape.Fmask, 
 *				mesh.Nv, mesh.x, mesh.y)
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
	int k,sk;
	for(k=0;k<K;k++){
		area [k] = 0.0;
		hmean[k] = 0.0;
		xmean[k] = 0.0;
		ymean[k] = 0.0;
		for(i=0;i<Np;i++){
			sk = (k*Np + i);
			hmean[k] += w[i]*J[sk]*h[sk];
			area [k] += w[i]*J[sk];
			xmean[k] += w[i]*J[sk]*x[sk];
			ymean[k] += w[i]*J[sk]*y[sk];
		}
		hmean[k] /= area[k];
		xmean[k] /= area[k];
		ymean[k] /= area[k];
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

	real alpha, a[4], b[2], gra[2];
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
			if (h[i] > hc*(1+TOTALERR) ){
				alpha = min(alpha, ( hvmax[v]-hc )/(h[i] - hc) );
			}else if(h[i] < hc*(1-TOTALERR)){
				alpha = min(alpha, ( hvmin[v]-hc )/(h[i] - hc) );
			}

			// mexPrintf("k=%d, v=%d, ind=%d, hi=%f, hmax=%f, hmin=%f, alpha=%f\n", 
				// k, f, v, h[i], hvmax[v], hvmin[v], alpha);
		}
		// mexPrintf("k=%d, alpha=%f\n",k, alpha);
		// get new veretx
		for(f=0;f<Nfaces;f++){
			i = k*Np + (int) Fmask[f]-1; // node index
			hv[f] = hc + alpha*(h[i] - hc);
		}
		// calculate the new gradient
		a[0] = xv[0] - xc; a[1] = yv[0] - yc;
		a[2] = xv[1] - xc; a[3] = yv[1] - yc;
		b[0] = hv[0] - hc; b[1] = hv[1] - hc;
		MatrixSolver2(a, b, gra);
		// mexPrintf("k=%d, phpx=%f, phpy=%f\n", k, gra[0], gra[1]);
		// new local variable
		GetLocalVar(Np, hc, xc, yc, x+k*Np, y+k*Np, gra[0], gra[1], hlim+k*Np);
	}

	free(w);
	free(hv); free(xv); free(yv);

	free(hmean); free(area);
	free(xmean); free(ymean);
	free(hvmax); free(hvmin);
}