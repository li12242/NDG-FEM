#include "Limiter.h"

/**
 * @brief
 * Use Hermite WENO scheme as a limiter to construct the limited result.
 * 
 * Usages:
 *		shape = mesh.Shape;
 * 		ind  = KXRCF_detector2d_Mex(h, u, v, mesh.J, mesh.sJ, 
 *							   shape.M, shape.Fmask, mesh.vmapM, mesh.vmapP, 
 *							   shape.Mes, mesh.nx, mesh.ny, distol)
 */

 void mexFunction (int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

 	/* check input & output */
	if (nrhs != 13)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	real *h    = mxGetPr(prhs[0]);
	real *u    = mxGetPr(prhs[1]);
	real *v    = mxGetPr(prhs[2]);
	real *J    = mxGetPr(prhs[3]);
	real *sJ   = mxGetPr(prhs[4]);
	real *M    = mxGetPr(prhs[5]);
	real *Fmask= mxGetPr(prhs[6]);
	real *vmapM= mxGetPr(prhs[7]);
	real *vmapP= mxGetPr(prhs[8]);
	real *Mes  = mxGetPr(prhs[9]);
	real *nx   = mxGetPr(prhs[10]);
	real *ny   = mxGetPr(prhs[11]);
	real *distol   = mxGetPr(prhs[12]);

	/* get dimensions */
	size_t Np, K;
	Np = mxGetM(prhs[0]); 
	K  = mxGetN(prhs[0]);
	size_t Nfaces,Nfp;
	Nfaces = mxGetM(prhs[6]);
	Nfp    = mxGetN(prhs[6]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)1, (mwSize)K, mxREAL);
	real *detector = mxGetPr(plhs[0]);

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
	// loop for all elements
	int k,f1,p;
	for(k=0;k<K;k++){
		real indicator = 0.0;
		real integral  = 0.0;
		real hmax;
		real radius, area = 0.0;

		hmax = fabs(h[k*Np]);
		for(p=0;p<Np;p++){
			int sk = (k*Np + p);
			hmax = max(hmax, fabs(h[sk]) );
			area += w[p]*J[sk];
		}
		radius = sqrt((5.0 - (real)Nfaces) *area);

		// mexPrintf("k=%d, hmax=%f, r=%f\n", k, hmax, radius);
		for(f1=0;f1<Nfaces;f1++){

            real flow_integral=0.0;
            for(p=0;p<Nfp;p++){
            	int sk = k*Np + (int) (*(Fmask+f1+p*Nfaces) - 1); // local index of face
            	int fk = k*Nfaces*Nfp+f1*Nfp+p;
            	real vn= nx[fk]*u[sk]+ny[fk]*v[sk];
            	flow_integral += ws[p]*sJ[fk]*vn;
            }

            // mexPrintf("k=%d, f=%d, Inflow=%f\n", k, f1, flow_integral);

            if(flow_integral < 0.0){

            	real face_len=0.0, sum=0.0, dQ;
            	for(p=0;p<Nfp;p++){
            		int fk  = k*Nfaces*Nfp+f1*Nfp+p;
            		int inM = (int) vmapM[fk] - 1;
            		int inP = (int) vmapP[fk] - 1;
            		
            		dQ   = h[inM] - h[inP];
            		sum += ws[p]*sJ[fk]*dQ;
            		face_len += ws[p]*sJ[fk];
            	}
            	indicator += sum;
            	integral  += face_len;
            }
		}

		// mexPrintf("k=%d, indicator=%f, integral=%f\n", k, indicator, integral);

		if(fabs(indicator) > (distol[0]*integral*hmax*radius) ){
			detector[k] = 1.0;
		}else{
			detector[k] = 0.0;
		}

		// mexPrintf("k=%d, detector=%f\n", k, detector[k]);
	}

	free(w); free(ws);
    return;
 }