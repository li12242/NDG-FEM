#include "Limiter.h"

#define TOTALERR 1e-12

real TVB_minmod(real a1, real a2, real dx, real factor);

/**
 * @brief
 * Use minmod function to limit the gradient and get the linear
 * limited result.
 *
 * Usages:
 *		shape = mesh.Shape;
 * 		hlim  = TVB2d_Mex...
 *			(h, mesh.J, mesh.sJ, shape.M, shape.Fmask, mesh.EToE, mesh.Mes, mesh.x, mesh.y, factor)
 */
void mexFunction (int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]){

	/* check input & output */
	if (nrhs != 10)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	real *h    = mxGetPr(prhs[0]);
	real *J    = mxGetPr(prhs[1]);
	real *sJ   = mxGetPr(prhs[2]);
	real *M    = mxGetPr(prhs[3]);
	real *Fmask= mxGetPr(prhs[4]);
	real *EToE = mxGetPr(prhs[5]);
	real *Mes  = mxGetPr(prhs[6]);
	real *x    = mxGetPr(prhs[7]);
	real *y    = mxGetPr(prhs[8]);
	real factor = mxGetScalar(prhs[9]);

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
	real *ws    = (real*) malloc(sizeof(real)*Nfp);

	/* volume/interface integral coefficient */
	int i,j,k,f,f1;
	for(i=0;i<Np;i++){
		w[i] = 0.0;
		for(j=0;j<Np;j++){
			w[i] += M[i*Np + j];
		}
	}

	for(i=0;i<Nfp;i++){
		ws[i] = 0.0;
		for(j=0;j<Np;j++){
			ws[i] += Mes[i*Np + j];
		}
	}
	// calculate volume mean value
	for(k=0;k<K;k++){
		real *fld[3], cmean[3];
		fld[0] = h+k*Np; 
		fld[1] = x+k*Np; 
		fld[2] = y+k*Np;
		cellMean(Np, 3, fld, w, J+k*Np, cmean, area+k);
		hmean[k] = cmean[0]; 
		xmean[k] = cmean[1]; 
		ymean[k] = cmean[2];
	}

	real *neigh_mean = (real*) malloc(sizeof(real)*Nfaces );
	real *neigh_xc   = (real*) malloc(sizeof(real)*Nfaces );
	real *neigh_yc   = (real*) malloc(sizeof(real)*Nfaces );
	real *xv         = (real*) malloc(sizeof(real)*Nfaces );
	real *yv         = (real*) malloc(sizeof(real)*Nfaces );
	// real *neigh_xf   = (real*) malloc(sizeof(real)*Nfaces );
	// real *neigh_yf   = (real*) malloc(sizeof(real)*Nfaces );
	// real *face_mean  = (real*) malloc(sizeof(real)*Nfaces );
	real *delta 	 = (real*) malloc(sizeof(real)*Nfaces );
	real alpha[2], a[4], df[2];

	for(k=0;k<K;k++){
		real xc = xmean[k];
		real yc = ymean[k];
		real hc = hmean[k];

		// get the adjacent cell averages
		real *fv[3], face_len;
		fv[0] = h+k*Np; 
		fv[1] = x+k*Np; 
		fv[2] = y+k*Np;
		for(f=0;f<Nfaces;f++){
			int e1 = (int) EToE[k+f*K] - 1;
			if(e1==k){
				real fmean[3];
				faceMean(Nfaces, Nfp, 3, fv, ws, sJ+k*Nfaces*Nfp+f*Nfp, Fmask+f, fmean, &face_len);
				neigh_mean[f] = fmean[0];
				neigh_xc[f]   = fmean[1];
				neigh_yc[f]   = fmean[2];
			}else{
				neigh_mean[f] = hmean[e1];
				neigh_xc[f] = xmean[e1];
				neigh_yc[f] = ymean[e1];
			}
		}

		// TVB limited vertex value
		for(f1=0;f1<Nfaces;f1++){
			int v  = k*Np + (int) (Fmask[f1]-1);
			int f2 = (int) (f1+Nfaces-1)%Nfaces;

			xv[f1] = x[v]; yv[f1] = y[v];
			df[0] = x[v] - xc; df[1] = y[v] - yc;
			a[0] = neigh_xc[f1] - xc; a[1] = neigh_xc[f2] - xc;
			a[2] = neigh_yc[f1] - yc; a[3] = neigh_yc[f2] - yc;

			matrixSolver2(a, df, alpha);

			real dhv = alpha[0]*(neigh_mean[f1]-hc)+alpha[1]*(neigh_mean[f2]-hc);
			real dx2 = sqrt(df[0]*df[0] + df[1]*df[1]);

			delta[f1] = TVB_minmod(h[v]-hc, dhv*1.5, dx2, factor);
			// mexPrintf("k=%d, v=%d, dh=%f, delta_h=%f, th=%f, delta=%f\n",
					// k, f1, h[v]-hc, dhv*v, dx2*factor[0], delta[f1]);
		}


	// for(k=0;k<K;k++){
	// 	real xc = xmean[k];
	// 	real yc = ymean[k];
	//
	// 	real *fv[3], fmean[3];
	// 	fv[0] = h+k*Np; fv[1] = x+k*Np; fv[2] = y+k*Np;
	// 	/* calculate average on face */
	// 	for(f=0;f<Nfaces;f++){
	// 		real face_len;
	// 		faceMean(Nfaces, Nfp, 3, fv, ws, sJ+k*Nfaces*Nfp+f*Nfp, Fmask+f, fmean, &face_len);
	// 		face_mean[f] = fmean[0]; neigh_xf[f]  = fmean[1]; neigh_yf[f]  = fmean[2];
	//
	// 		// mexPrintf("k=%d, Nfaces=%d, xf=%f, yf=%f, face_len=%f, face_mean=%f\n",
	// 		// 		k, f, neigh_xf[f], neigh_yf[f], face_len, face_mean[f]);
	// 		// averages of adjacent elements
	// 		int e1 = (int)EToE[k+K*f] - 1;
	// 		if(e1 == k){ // for boundary
	// 			neigh_mean[f] = face_mean[f];
	// 			neigh_xc[f] = neigh_xf[f];
	// 			neigh_yc[f] = neigh_yf[f];
	// 		}else{
	// 			neigh_mean[f] = hmean[e1];
	// 			neigh_xc[f] = xmean[e1];
	// 			neigh_yc[f] = ymean[e1];
	// 		}
	// 	}
	//
	// 	for(f1=0;f1<Nfaces;f1++){
	// 		/* find best other face */
	// 		real max_alpha = -1.0;
	// 		int  best_face = -1;
	// 		real beat_alpha[2];
	//
	// 		a[0]  = neigh_xc[f1] - xc;
	// 		a[2]  = neigh_yc[f1] - yc;
	// 		df[0] = neigh_xf[f1] - xc;
	// 		df[1] = neigh_yf[f1] - yc;
	//
	// 		for(f2=0;f2<Nfaces;f2++){
	// 			if(f1==f2)
	// 				continue;
	//
	// 			/* calculate of alpha */
	// 			a[1] = neigh_xc[f2] - xc;
	// 			a[3] = neigh_yc[f2] - yc;
	//
	// 			matrixSolver2(a, df, alpha);
	//
	// 			real alpha_det = sqrt(alpha[0]*alpha[0] + alpha[1]*alpha[1]);
	// 			if ((alpha[0] > -TOTALERR) & (alpha[1] > -TOTALERR) & alpha[0]/alpha_det > max_alpha ){
	// 				best_face = f2;
	// 				max_alpha = alpha[0]/alpha_det;
	// 				beat_alpha[0] = alpha[0];
	// 				beat_alpha[1] = alpha[1];
	// 			}
	// 		}
	//
	// 		real dhc, dhf, len;
	//
	// 		dhc    = alpha[0]*(neigh_mean[f1] - hmean[k]) + alpha[1]*(neigh_mean[best_face] - hmean[k]);
	// 		dhf    = face_mean[f1] - hmean[k];
	// 		len    = sqrt( (neigh_xc[f1]-xc)*(neigh_xc[f1]-xc)+(neigh_yc[f1]-yc)*(neigh_yc[f1]-yc) );
	//
	// 		delta[f1] = TVB_minmod(dhf, dhc, len, *factor);
	//
	// 		// mexPrintf("k=%d, Nfaces=%d, f1=%d, xf=%f, yf=%f, delta=%f\n",
	// 				// k, f1, best_face, neigh_xf[f1], neigh_yf[f1], delta[f1]);
	// 	}
	// 	/* correct the average to sum to 0.0 */
		real sum = 0.0;
		for(f=0;f<Nfaces;f++){
			sum += delta[f];
		}

		if( abs(sum)>TOTALERR ){
			real pos = 0.0;
			real neg = 0.0;

			for(f=0;f<Nfaces;f++){
				pos += max( delta[f], 0.0);
				neg += max(-delta[f], 0.0);
			}

			for(f=0;f<Nfaces;f++){
				delta[f] = min(1.0, neg/pos)*max(delta[f], 0.0)
					- min(1.0, pos/neg)*max(-delta[f], 0.0);
			}
		}

		for(f=0;f<Nfaces;f++){
			delta[f] += hmean[k];
		}

		/* reconstruct the cell value */
		real qpx = 0.0;
        real qpy = 0.0;

		VA_meanGradient(Nfaces, xv, yv, delta, xc, yc, hc, &qpx, &qpy);
        getLocalVar(Np, hc, xc, yc, x+k*Np, y+k*Np, qpx, qpy, hlim+k*Np);

	}

	free(delta);
	free(neigh_xc);
	free(neigh_yc);
	free(xv);
	free(yv);
	// free(neigh_xf);
	// free(neigh_yf);
	// free(face_mean);
	free(neigh_mean);

	free(hmean);
	free(area);
	free(xmean);
	free(ymean);
	free(w);
	free(ws);
	return;
}

/*
 * modified minmod function of TVB limiter.
 * m = TVB_minmod(a1, a2)
 */
real TVB_minmod(real a1, real a2, real dx2, real TVB_factor){
	real m;
	// mexPrintf("a1=%f, a2=%f, Th=%f\n", a1, a2, TVB_factor*dx2);
	if(fabs(a1)<TVB_factor*dx2){
		m = a1;
		// mexPrintf("a1=%d, Th=%f\n", fabs(a1), TVB_factor*dx2);
	}else{
		real a[2];
		a[0] = a1; a[1] = a2;
		minmod(2, a, &m);
		// mexPrintf("a1=%f, a2=%f, minmod=%f\n", a[0], a[1], m);
	}
	return m;
}
