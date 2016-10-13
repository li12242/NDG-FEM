#include "Limiter.h"

#define TOTALERR 1e-12

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

	if(Nfaces!=4){
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
	int i,j,k,f;
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

	// double alpha[2], a[4], b[2]; // coefficient for equations

	for(k=0;k<K;k++){
		double xc = xmean[k];
		double yc = ymean[k];
		double hc = hmean[k];
		double face_len;
		double xf[4], yf[4], hf[4]; // value of middle point on face
		double he[4]; // value of adjacent element
		double delta[4]; // limited middle point value
	
		double *p[3], fmean[3];
		p[0] = h+k*Np; 
		p[1] = x+k*Np; 
		p[2] = y+k*Np;
		
		for(f=0;f<Nfaces;f++){
			/* calculate average on face */
			faceMean(Nfaces, Nfp, 3, p, 
				ws, sJ+k*Nfaces*Nfp+f*Nfp, 
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
				// xe[f] = 2*xf[f] - xc;
				// ye[f] = 2*yf[f] - yc;
			}else{
				he[f] = hmean[e1];
				// xe[f] = xmean[e1];
				// ye[f] = ymean[e1];
			}
		}

		// double sum = 0.0;
		// for(f=0;f<Nfaces;f++){
		// 	sum += hf[f];
		// }

		// if(fabs(sum - 4.0*hc) > -TOTALERR)
		// 	mexPrintf("Warning, k=%d the average value of face mean is not equal to hc\n", k);

		/* limit the face middle point */
		for(f=0;f<Nfaces;f++){
			int e1 = f;
			int e2 = (f+2)%Nfaces;

			double dh[3];
			dh[0] = hf[f] - hc;
			dh[1] = (he[e1] - hc)*factor;
			dh[2] = (hc - he[e2])*factor;

			minmod(3, dh, delta+f);
		}

		/* reconstruct the cell value */
		double qpx = 0.0;
        double qpy = 0.0;

        VA_meanGradient(Nfaces, xf, yf, delta, xc, yc, hc, &qpx, &qpy);
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