#include "Limiter.h"

/**
 * @brief
 * Use Hermite WENO scheme as a limiter to construct the limited result.
 * 
 * Usages:
 *		shape = mesh.Shape;
 * 		hlim  = TVB2d_Mex...
 *			(h, mesh.J, mesh.sJ, shape.M, shape.Fmask, mesh.EToE, 
 *				mesh.Mes, mesh.x, mesh.y, eps_o, eps_w, gam0)
 */

 void mexFunction (int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

 	/* check input & output */
	if (nrhs != 12)
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
	real *eps_o    = mxGetPr(prhs[9]);
	real *eps_w    = mxGetPr(prhs[10]);
	real *gam0     = mxGetPr(prhs[11]);

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

	int k,sk,sp,f1,f2,e1,e2,fnp;
	/* cell averages */
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

		/* gradient */
		phpx[k] = 0.0;
		phpy[k] = 0.0;
		for(f1=0;f1<Nfaces;f1++){
			/* vertex index of face */
            int l1 = k*Np + (int) (Fmask[f1 + 0*Nfaces] - 1);
            int l2 = k*Np + (int) (Fmask[f1 + (Nfp-1)*Nfaces] - 1);
            
            /* mean value on edge */
            real dx = (real) (x[l2] - x[l1]);
            real dy = (real) (y[l2] - y[l1]);
            real face_len, face_mean;
			FaceMean(Nfaces, Nfp, h+k*Np, ws, sJ+k*Nfaces*Nfp+f1*Nfp, Fmask+f1,
				&face_mean, &face_len);

			/* Green-Gauss formulation */
			phpx[k] += face_mean*dy;
			phpy[k] -= face_mean*dx;
		}
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
			/* get adjacent info */
			if(e1 == k){ // use face info
				real face_mean;
				FaceMean(Nfaces, Nfp, h+k*Np, ws, sJ+k*Nfaces*Nfp+f1*Nfp, Fmask+f1,
					&h1, &face_len);
				FaceMean(Nfaces, Nfp, x+k*Np, ws, sJ+k*Nfaces*Nfp+f1*Nfp, Fmask+f1,
					&x1, &face_len);
				FaceMean(Nfaces, Nfp, y+k*Np, ws, sJ+k*Nfaces*Nfp+f1*Nfp, Fmask+f1,
					&y1, &face_len);
			}else{ // use e1
				x1 = xmean[e1];
				y1 = ymean[e1];
				h1 = hmean[e1];
			}

			if(e2 ==k ){ // use face info
				real face_mean;
				FaceMean(Nfaces, Nfp, h+k*Np, ws, sJ+k*Nfaces*Nfp+f2*Nfp, Fmask+f2,
					&h2, &face_len);
				FaceMean(Nfaces, Nfp, x+k*Np, ws, sJ+k*Nfaces*Nfp+f2*Nfp, Fmask+f2,
					&x2, &face_len);
				FaceMean(Nfaces, Nfp, y+k*Np, ws, sJ+k*Nfaces*Nfp+f2*Nfp, Fmask+f2,
					&y2, &face_len);
			}else{
				x2 = xmean[e2];
				y2 = ymean[e2];
				h2 = hmean[e2];
			}

			a[0] = x1 - xc; a[1] = y1 - yc;
			a[2] = x2 - xc; a[3] = y2 - yc;
			f[0] = h1 - hc; f[1] = h2 - hc;
			MatrixSolver2(a, f, grd);
			phpx_stencil[f1*2+1] = grd[0];
			phpy_stencil[f1*2+1] = grd[1];

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

				MatrixSolver2(a, f, grd);
				phpx_stencil[f1*2+2] = grd[0];
				phpy_stencil[f1*2+2] = grd[1];
			}else{ //use adjacent cell gradient
				phpx_stencil[f1*2+2] = phpx[e1];
				phpy_stencil[f1*2+2] = phpy[e1];
			}
		}

		real sum_w = 0.0;
		/* oscillation indicator and weights */
		for(i=0;i<(2*Nfaces+1);i++){
			polys_o[i] = phpx_stencil[i]*phpx_stencil[i] 
					+ phpy_stencil[i]*phpy_stencil[i];
			real sum = 0;
			real *hlocal = (real*) malloc(sizeof(real)*Np );
			GetLocalVar(Np, hc, xc, yc, x+k*Np, y+k*Np, 
				phpx_stencil[i], phpy_stencil[i], hlocal);
			for(j=0;j<Np;j++){
				sk = (k*Np + i);
				sum += w[i]*J[sk]*(hlocal[j]*hlocal[j]+eps_o[0]);
			}

			polys_o[i] /= sum;
			polys_w[i] = 1.0*pow(eps_w[0] + polys_o[i], -gam0[0]);
			sum_w += polys_w[i];
		}

		real phpxlim, phpylim;
		for(i=0;i<(2*Nfaces+1);i++){
			polys_w[i] /= sum_w;

			phpxlim += polys_w[i]*phpx_stencil[i];
			phpylim += polys_w[i]*phpy_stencil[i];
		}

		GetLocalVar(Np, hc, xc, yc, x+k*Np, y+k*Np, 
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