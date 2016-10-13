#include "Limiter.h"

/**
 * @brief
 * Use Hermite WENO scheme as a limiter to construct the limited result.
 *
 * Usages:
 *		shape = mesh.Shape;
 * 		hlim  = HWENO2d_Mex(h, mesh.J, mesh.sJ, shape.M, shape.Fmask, mesh.EToE,
 *							   mesh.Mes, mesh.x, mesh.y, eps_w, gam0)
 */

 void mexFunction (int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]){

 	/* check input & output */
	if (nrhs != 11)
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
	real *eps_w    = mxGetPr(prhs[9]);
	real *gam0     = mxGetPr(prhs[10]);

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

	//elemental integral coefficient
	real *w     = (real*) malloc(sizeof(real)*Np);
	real *ws    = (real*) malloc(sizeof(real)*Nfp);

	/* volume/interface integral coefficient */
	int i,j;
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

	/* cell averages */
	real *hmean = (real*) malloc(sizeof(real)*K );
	real *area  = (real*) malloc(sizeof(real)*K );
	real *xmean = (real*) malloc(sizeof(real)*K );
	real *ymean = (real*) malloc(sizeof(real)*K );
	real *phpx  = (real*) malloc(sizeof(real)*K );
	real *phpy  = (real*) malloc(sizeof(real)*K );

	real *phpx_stencil = (real*) malloc(sizeof(real)*(2*Nfaces+1) );
	real *phpy_stencil = (real*) malloc(sizeof(real)*(2*Nfaces+1) );
	real *polys_o = (real*) malloc(sizeof(real)*(2*Nfaces+1));
	real *polys_w = (real*) malloc(sizeof(real)*(2*Nfaces+1));

	int k,f1,f2,e1,e2;
	/* cell averages */
	for(k=0;k<K;k++){
        real *fld[3], cmean[3];
		fld[0] = h+k*Np; fld[1] = x+k*Np; fld[2] = y+k*Np;
		cellMean(Np, 3, fld, w, J+k*Np, cmean, area+k);
		hmean[k] = cmean[0];
		xmean[k] = cmean[1];
		ymean[k] = cmean[2];

		/* gradient */
		phpx[k] = 0.0;
		phpy[k] = 0.0;
        real *hv[1];
        hv[0] = h+k*Np;
		for(f1=0;f1<Nfaces;f1++){
			/* vertex index of face */
            int l1 = k*Np + (int) (Fmask[f1 + 0*Nfaces] - 1);
            int l2 = k*Np + (int) (Fmask[f1 + (Nfp-1)*Nfaces] - 1);

            /* mean value on edge */
            real dx = (real) (x[l2] - x[l1]);
            real dy = (real) (y[l2] - y[l1]);
            real face_len, face_mean;

 			faceMean(Nfaces, Nfp, 1, hv, ws, sJ+k*Nfaces*Nfp+f1*Nfp, Fmask+f1,
				&face_mean, &face_len);

			/* Green-Gauss formulation */
			phpx[k] += face_mean*dy;
			phpy[k] -= face_mean*dx;

			// mexPrintf("k=%d, f=%d, f_mean=%f\n", k, f1, face_mean);
		}
		// mexPrintf("k=%d, phpx=%f, phpy=%f\n", k, phpx[k], phpy[k]);
	}

	for(k=0;k<K;k++){
		phpx_stencil[0] = phpx[k];
		phpy_stencil[0] = phpy[k];
		real xc = xmean[k];
		real yc = ymean[k];
		real hc = hmean[k];

		for(f1=0;f1<Nfaces;f1++){
			f2 = (f1+1)%Nfaces;
			e1 = (int) (int)EToE[k+f1*K] - 1;
			e2 = (int) (int)EToE[k+f2*K] - 1;

            real x1, x2, y1, y2, h1, h2, face_len;
			real a[4], f[2], grd[2];
            real *fv[3],fmean[3];
            fv[0] = h+k*Np; fv[1] = x+k*Np; fv[2] = y+k*Np;
			/* get adjacent info */
			if(e1 == k){ // use face info

				faceMean(Nfaces, Nfp, 3, fv, ws, sJ+k*Nfaces*Nfp+f1*Nfp, Fmask+f1,
					fmean, &face_len);
                h1 = fmean[0]; x1 = fmean[1]; y1 = fmean[2];
				// FaceMean(Nfaces, Nfp, x+k*Np, ws, sJ+k*Nfaces*Nfp+f1*Nfp, Fmask+f1,
				// 	&x1, &face_len);
				// FaceMean(Nfaces, Nfp, y+k*Np, ws, sJ+k*Nfaces*Nfp+f1*Nfp, Fmask+f1,
				// 	&y1, &face_len);
			}else{ // use e1
				x1 = xmean[e1];
				y1 = ymean[e1];
				h1 = hmean[e1];
			}

			if(e2 == k ){ // use face info
                faceMean(Nfaces, Nfp, 3, fv, ws, sJ+k*Nfaces*Nfp+f2*Nfp, Fmask+f2,
					fmean, &face_len);
                h2 = fmean[0]; x2 = fmean[1]; y2 = fmean[2];

				// FaceMean(Nfaces, Nfp, h+k*Np, ws, sJ+k*Nfaces*Nfp+f2*Nfp, Fmask+f2,
				// 	&h2, &face_len);
				// FaceMean(Nfaces, Nfp, x+k*Np, ws, sJ+k*Nfaces*Nfp+f2*Nfp, Fmask+f2,
				// 	&x2, &face_len);
				// FaceMean(Nfaces, Nfp, y+k*Np, ws, sJ+k*Nfaces*Nfp+f2*Nfp, Fmask+f2,
				// 	&y2, &face_len);
			}else{
				x2 = xmean[e2];
				y2 = ymean[e2];
				h2 = hmean[e2];
			}

			a[0] = x1 - xc; a[1] = y1 - yc;
			a[2] = x2 - xc; a[3] = y2 - yc;
			f[0] = h1 - hc; f[1] = h2 - hc;
			matrixSolver2(a, f, grd);
			phpx_stencil[f1*2+1] = grd[0];
			phpy_stencil[f1*2+1] = grd[1];

			// mexPrintf("k=%d, phpx[%d]=%f, phpy[%d]=%f\n", k, f1*2+1, phpx_stencil[f1*2+1], f1*2+1, phpy_stencil[f1*2+1]);

			if (e1 == k){ //use face vertex and centre to calculate the gradient
				/* vertex index of face */
	            int l1 = k*Np + (int) (Fmask[f1 + 0*Nfaces] - 1);
	            int l2 = k*Np + (int) (Fmask[f1 + (Nfp-1)*Nfaces] - 1);

	            x1 = (real) x[l1]; x2 = (real) x[l2];
	            y1 = (real) y[l1]; y2 = (real) y[l2];
	            h1 = (real) h[l1]; h2 = (real) h[l2];

	            a[0] = x1 - xc; a[1] = y1 - yc;
				a[2] = x2 - xc; a[3] = y2 - yc;
				f[0] = h1 - hc; f[1] = h2 - hc;

				matrixSolver2(a, f, grd);
				phpx_stencil[f1*2+2] = grd[0];
				phpy_stencil[f1*2+2] = grd[1];
			}else{ //use adjacent cell gradient
				phpx_stencil[f1*2+2] = phpx[e1];
				phpy_stencil[f1*2+2] = phpy[e1];
			}

			// mexPrintf("k=%d, phpx[%d]=%f, phpy[%d]=%f\n", k, f1*2+2, phpx_stencil[f1*2+2], f1*2+2, phpy_stencil[f1*2+2]);
		}

		real sum_w = 0.0;
		real *hlocal = (real*) malloc(sizeof(real)*Np );
		/* oscillation indicator and weights */
		for(i=0;i<(2*Nfaces+1);i++){
			polys_o[i] = phpx_stencil[i]*phpx_stencil[i]
					+ phpy_stencil[i]*phpy_stencil[i];

			getLocalVar(Np, hc, xc, yc, x+k*Np, y+k*Np,
				phpx_stencil[i], phpy_stencil[i], hlocal);

			// mexPrintf("hlocal=%f, %f, %f, %f\n", hlocal[0], hlocal[1], hlocal[2], hlocal[3] );
			// real sum = 0;
			// for(j=0;j<Np;j++){
			// 	sk = (k*Np + j);
			// 	sum += w[j]*J[sk]*(hlocal[j]*hlocal[j]+eps_o[0]);
			// }
			// sum /= area[k];

			polys_o[i] /= area[k];
			polys_o[i] = sqrt(polys_o[i]);
			polys_w[i] = 1.0*pow(eps_w[0] + polys_o[i], -gam0[0]);
			sum_w += polys_w[i];

			// mexPrintf("k=%d, sum=%f, polys_o[%d]=%f\n", k, sum, i, polys_o[i]);
		}

		free(hlocal);
		real phpxlim=0.0, phpylim=0.0;
		for(i=0;i<(2*Nfaces+1);i++){
			polys_w[i] /= sum_w;

			phpxlim += polys_w[i]*phpx_stencil[i];
			phpylim += polys_w[i]*phpy_stencil[i];

			// mexPrintf("k=%d, polys_w[%d]=%f\n", k, i, polys_w[i]);
		}

		// mexPrintf("k=%d, phpxlim=%f, phpylim=%f\n", k, phpxlim, phpylim);

		getLocalVar(Np, hc, xc, yc, x+k*Np, y+k*Np,
			phpxlim, phpylim, hlim+k*Np);

	}

	free(w);
	free(ws);

	free(hmean);
	free(area);
	free(xmean);
	free(ymean);
	free(phpx);
	free(phpy);
	free(phpx_stencil);
	free(phpy_stencil);

	free(polys_o);
	free(polys_w);

	return;

 }
