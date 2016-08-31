#include "Limiter.h"


/**
 * @brief
 * Use minmod function to limit the gradient and get the linear 
 * limited result.
 * 
 * Usages:
 *		shape = mesh.Shape;
 * 		flag  = TVB_detector2d_Mex(...
 *						h, mesh.J, shape.M, shape.Fmask, mesh.EToE, factor)
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
	real *factor = mxGetPr(prhs[5]);

	/* get dimensions */
	size_t Np, K;
	Np = mxGetM(prhs[0]); 
	K  = mxGetN(prhs[0]);
	size_t Nfaces,Nfp;
	Nfaces = mxGetM(prhs[3]);
	Nfp    = mxGetN(prhs[3]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)1, (mwSize)K, mxREAL);
	real *detector = mxGetPr(plhs[0]);

	/* cell averages */
	real *hmean = (real*) malloc(sizeof(real)*K );
	real *area  = (real*) malloc(sizeof(real)*K );
	//elemental integral coefficient
	real *w     = (real*) malloc(sizeof(real)*Np); 

	/* volume/interface integral coefficient */
	int i,j;
	for(i=0;i<Np;i++){
		w[i] = 0.0;
		for(j=0;j<Np;j++){
			w[i] += M[i*Np + j];
		}
	}

	// calculate volume mean value
	int sk,k,f,p;
	for(k=0;k<K;k++){
		area [k] = 0.0;
		hmean[k] = 0.0;
		for(i=0;i<Np;i++){
			sk = (k*Np + i);
			hmean[k] += w[i]*J[sk]*h[sk];
			area [k] += w[i]*J[sk];
		}
		hmean[k] /= area[k];
	}

	for(k=0;k<K;k++){
		detector[k] = 0.0;
		
		if(detector[k]>1.0e-10)
			continue;

		for(f=0;f<Nfaces;f++){
			int e2 = (int) EToE[k+f*K]-1;
			if(e2 == k)
				continue;

			real face_max = h[k*Np + (int)(*(Fmask+f) - 1)];
			real face_min = h[k*Np + (int)(*(Fmask+f) - 1)];
			for(p=0;p<Nfp;p++){
				int sk = k*Np + (int) (*(Fmask+f+p*Nfaces) - 1);
				face_max = max(face_max, h[sk]);
				face_min = min(face_min, h[sk]);
			}

			if(face_max > max(hmean[k], hmean[e2])+factor[0] ) 
				detector[k] = 1.0;
			if(face_min > min(hmean[k], hmean[e2])+factor[0] ) 
				detector[k] = 1.0;

			// mexPrintf("k=%d, f=%d, fmax=%f, fmin=%f, e1=%f, e2=%f, ind=%f\n", 
				// k, f, face_max, face_min, hmean[k], hmean[e2], detector[k]);
		}
	}

	free(hmean); free(area); free(w);
	return;
}