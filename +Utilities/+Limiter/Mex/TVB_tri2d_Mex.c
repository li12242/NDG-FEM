#include "Limiter.h"

#define TOTALERR 1e-12

double TVB_minmod(double a1, double a2, double dx, double factor);

/**
 * @brief
 * Use minmod function to limit the gradient and get the linear
 * limited result.
 *
 * Usages:
 *		shape = mesh.Shape;
 * 		hlim  = TVB2d_Mex...
 *			(h, mesh.J, mesh.sJ, shape.M, shape.Fmask, 
 *				mesh.EToE, mesh.Mes, mesh.x, mesh.y, factor)
 */
void mexFunction (int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]){

	/* check input & output */
	if (nrhs != 10)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	double *h    = mxGetPr(prhs[0]);
	double *J    = mxGetPr(prhs[1]);
	double *sJ   = mxGetPr(prhs[2]);
	double *M    = mxGetPr(prhs[3]);
	double *Fmask= mxGetPr(prhs[4]);
	double *EToE = mxGetPr(prhs[5]);
	double *Mes  = mxGetPr(prhs[6]);
	double *x    = mxGetPr(prhs[7]);
	double *y    = mxGetPr(prhs[8]);
	double factor = mxGetScalar(prhs[9]);

	/* get dimensions */
	size_t Np, K;
	Np = mxGetM(prhs[0]);
	K  = mxGetN(prhs[0]);
	size_t Nfaces,Nfp;
	Nfaces = mxGetM(prhs[4]);
	Nfp    = mxGetN(prhs[4]);

	if(Nfaces!=3){
		mexPrintf("Wrong number of Nfaces=%d for TVB_tri limiter", Nfaces);
		exit(-1);
	}

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	double *hlim = mxGetPr(plhs[0]);

	/* cell averages */
	double *hmean = (double*) malloc(sizeof(double)*K );
	double *area  = (double*) malloc(sizeof(double)*K );
	double *xmean = (double*) malloc(sizeof(double)*K );
	double *ymean = (double*) malloc(sizeof(double)*K );
	//elemental integral coefficient
	double *w     = (double*) malloc(sizeof(double)*Np);
	double *ws    = (double*) malloc(sizeof(double)*Nfp);

	/* volume/interface integral coefficient */
	int i,j,k,f,f1,f2;
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
	/* calculate volume mean value */
	for(k=0;k<K;k++){
		double *fld[3], cmean[3];
		fld[0] = h+k*Np; 
		fld[1] = x+k*Np; 
		fld[2] = y+k*Np;
		cellMean(Np, 3, fld, w, J+k*Np, cmean, area+k);
		hmean[k] = cmean[0]; 
		xmean[k] = cmean[1]; 
		ymean[k] = cmean[2];
	}
    
	double alpha[2], a[4], b[2];

	for(k=0;k<K;k++){
		double xc = xmean[k];
		double yc = ymean[k];
		double hc = hmean[k];
		double face_len;
		double xf[3], yf[3], hf[3]; // value of middle point on face
		double xe[3], ye[3], he[3]; // value of adjacent element
		double delta[3]; // limited middle point value
	
		double *p[3], fmean[3];
		p[0] = h+k*Np; 
		p[1] = x+k*Np; 
		p[2] = y+k*Np;
		
		for(f=0;f<Nfaces;f++){
			/* calculate average on face */
			faceMean(Nfaces, Nfp, 3, p, ws, sJ+k*Nfaces*Nfp+f*Nfp, 
				Fmask+f, fmean, &face_len);
			hf[f] = fmean[0]; 
			xf[f] = fmean[1]; 
			yf[f] = fmean[2];

			/* value of adjacent elements */
			int e1 = (int)EToE[k+K*f] - 1;
			if(e1 == k){ // for boundary
				// he[f] = hf[f];
				// xe[f] = xf[f];
				// ye[f] = yf[f];
				he[f] = hmean[e1];
				xe[f] = 2*xf[f] - xc;
				ye[f] = 2*yf[f] - yc;
			}else{
				he[f] = hmean[e1];
				xe[f] = xmean[e1];
				ye[f] = ymean[e1];
			}
		}
	
		
		for(f1=0;f1<Nfaces;f1++){
			/* find best other face */
			double max_alpha = -1.0;
			int    best_face = -1;
			double best_alpha[2];
	
			a[0] = xe[f1] - xc;
			a[2] = ye[f1] - yc;
			b[0] = xf[f1] - xc;
			b[1] = yf[f1] - yc;
	
			for(f2=0;f2<Nfaces;f2++){
				if(f1==f2)
					continue;
	
				/* calculate of alpha */
				a[1] = xe[f2] - xc;
				a[3] = ye[f2] - yc;
	
				matrixSolver2(a, b, alpha);
	
				double alpha_det = sqrt(alpha[0]*alpha[0] + alpha[1]*alpha[1]);
				if ((alpha[0] > -TOTALERR) & (alpha[1] > -TOTALERR) 
							& alpha[0]/alpha_det > max_alpha ){

					best_face = f2;
					max_alpha = alpha[0]/alpha_det;
					best_alpha[0] = alpha[0];
					best_alpha[1] = alpha[1];
				}
			}
			// mexPrintf("k=%d, f=%d, f2=%d, alpha=[%f, %f]\n", 
			// 	k,f1,best_face,best_alpha[0],best_alpha[1]);
	
			double dhe, dhf, dx2;
	
			dhe    = best_alpha[0]*(he[f1] - hc) 
						+ best_alpha[1]*(he[best_face] - hc);
			dhf    = hf[f1] - hc;
			dx2    = (xe[f1]-xc)*(xe[f1]-xc)+(ye[f1]-yc)*(ye[f1]-yc);
	
			delta[f1] = TVB_minmod(dhf, 1.5*dhe, dx2, factor);

			// mexPrintf("k=%d, f=%d, du=%f, due=%f, M=%f, delta=%f\n",
			// 	k,f1,dhf,dhe,factor*dx2,delta[f1]);
		}

		/* correct the average to sum to 0.0 */
		double sum = 0.0;
		for(f=0;f<Nfaces;f++){
			sum += delta[f];
		}

		if( abs(sum)>TOTALERR ){
			double pos = 0.0;
			double neg = 0.0;

			for(f=0;f<Nfaces;f++){
				pos += max( delta[f], 0.0);
				neg += max(-delta[f], 0.0);
			}

			for(f=0;f<Nfaces;f++){
				delta[f] = min(1.0, neg/pos)*max(delta[f], 0.0)
					- min(1.0, pos/neg)*max(-delta[f], 0.0);
			}
			// mexPrintf("k=%d, pos=%e, neg=%e, delta=%f\n", k, pos, neg, delta);
		}

		for(f=0;f<Nfaces;f++){
			delta[f] += hmean[k];
		}

		/* reconstruct the cell value */
		double qpx = 0.0;
        double qpy = 0.0;

		a[0] = xf[0] - xf[1]; a[1] = yf[0] - yf[1]; b[0] = delta[0] - delta[1];
		a[2] = xf[0] - xf[2]; a[3] = yf[0] - yf[2]; b[1] = delta[0] - delta[2];

		matrixSolver2(a, b, alpha);
		
		qpx = alpha[0];
		qpy = alpha[1];

        getLocalVar(Np, hc, xc, yc, x+k*Np, y+k*Np, qpx, qpy, hlim+k*Np);

	}

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
double TVB_minmod(double a1, double a2, double dx2, double TVB_factor){
	double m;
	// mexPrintf("a1=%f, a2=%f, Th=%f\n", a1, a2, TVB_factor*dx2);
	if(fabs(a1)<TVB_factor*dx2){
		m = a1;
		// mexPrintf("a1=%d, Th=%f\n", fabs(a1), TVB_factor*dx2);
	}else{
		double a[2];
		a[0] = a1; a[1] = a2;
		minmod(2, a, &m);
		// mexPrintf("a1=%f, a2=%f, minmod=%f\n", a[0], a[1], m);
	}
	return m;
}
