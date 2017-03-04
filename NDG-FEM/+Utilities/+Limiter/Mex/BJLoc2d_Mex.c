#include "mex.h"

#define real double

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

/*
 *
 * Usages:
 * 		hlim = SLLoc2d_Mex...
 *         (h, mesh.J, mesh.Shape.M, mesh.Shape.Fmask, mesh.EToE, mesh.x, mesh.y, beta)
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
	real *Fmask= mxGetPr(prhs[3]);
	real *EToE = mxGetPr(prhs[4]);
	real *x    = mxGetPr(prhs[5]);
	real *y    = mxGetPr(prhs[6]);
	real beta  = mxGetScalar(prhs[7]);

	/* get dimensions */
	size_t Np, K;
	Np = mxGetM(prhs[0]); 
	K  = mxGetN(prhs[0]);
	size_t Nfaces,Nfp;
	Nfp    = mxGetM(prhs[3]);
	Nfaces = mxGetN(prhs[3]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	real *hlim = mxGetPr(plhs[0]);

	/* get cell average */
	real *area  = (real*) malloc(sizeof(real)*K );
	real *hmean = (real*) malloc(sizeof(real)*K );

	//elemental integral coefficient
	real *w     = (real*) malloc(sizeof(real)*Np); 
	int i,j,k,sk;
	for(i=0;i<Np;i++){
		w[i] = 0.0;
		for(j=0;j<Np;j++){
			w[i] += M[i*Np + j];
		}
	}
	// cell area and scalar averages
	sk = 0;
	for (i=0;i<K;i++){
		area[i]  = 0.0;
		hmean[i] = 0.0;
		for(j=0;j<Np;j++){
			hmean[i] += w[j]*J[sk]*h[sk];
			area[i]  += w[j]*J[sk];
			sk++;
		}
		hmean[i] /= area[i];
	}

	/* get unlimited gradient */
	real phpx, phpy; // gradient
	for (k=0;k<K;k++){
		phpx = 0.0;
		phpy = 0.0;

		real hmin = hmean[k];
		real hmax = hmean[k];

		for(i=0;i<Nfaces;i++){
			int e = (int)EToE[k+i*K] - 1;

			hmin = min(hmin, hmean[e]);
			hmax = max(hmax, hmean[e]);

			// int v1 = (int) Fmask[i] - 1;
			// int v2 = (int) Fmask[i+(Nfp-1)*Nfaces] - 1;
			int v1 = (int) Fmask[Nfp*i] - 1; // local vertex indics
			int v2 = (int) Fmask[Nfp*i + (Nfp-1)] - 1;

			int sk1 = k*Np + v1;
			int sk2 = k*Np + v2;
			real hmean1 = (h[sk1] + h[sk2])/2.0;

			real dx = x[sk2] - x[sk1];
			real dy = y[sk2] - y[sk1];

			phpx += hmean1*dy;
			phpy -= hmean1*dx;
		}

		phpx /= area[k];
		phpy /= area[k];

		real xc=0.0,yc=0.0;
        real t, psi = 1.0;
        for(i=0;i<Np;i++){
            sk = (k*Np + i);
            real hval = h[sk], rk = 1.0;
            /* compute the limiter `psi` */
            if(hval > hmean[k]){
                rk = (hmax - hmean[k])/(hval - hmean[k]);
            }else if(hval < hmean[k]){
                rk = (hmin - hmean[k])/(hval - hmean[k]);
            }
            t   = max(min(beta*rk, 1.0), min(rk, beta));
            psi = min(psi, t);

            /* compute the centre of cell */
            xc += w[i]*J[sk]*x[sk];
            yc += w[i]*J[sk]*y[sk];
        }
        xc /= area[k];
        yc /= area[k];

        /* compute the limited gradient */
        phpx *= psi;
        phpy *= psi;

        /* reconstruction of each element */
        for(i=0;i<Np;i++){
            sk = (k*Np + i);
            real dx = (real)(x[sk] - xc);
            real dy = (real)(y[sk] - yc);

            hlim[sk] = hmean[k] + dx*phpx + dy*phpy;
        }

	}

	free(w);
	free(area);
	free(hmean);
	return;
}