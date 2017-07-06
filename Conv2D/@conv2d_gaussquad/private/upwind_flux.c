#include "conv2d.h"

void upwind_flux(double f_M, double f_P, double uM, double vM, 
	double nx, double ny, double *numflux){
	
	const double unM = uM*nx + vM*ny;
	if (unM > 0){
		*numflux = f_M*unM;
	}else{
		*numflux = f_P*unM;
	}
	return;
}

/* @brief calculate the surface flux deviation for strong form.
 * 
 * Usages:
 * 	[dflux] = upwind_flux(h, h_ext, u, v, nx, ny, eidM, eidP, eidtype);
 */
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]){

	/* check input & output */
	if (nrhs != 9) mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1) mexErrMsgTxt("Wrong number of output arguments.");

	/* get inputs */
	double *h = mxGetPr(prhs[0]);
    double *h_ext = mxGetPr(prhs[1]);
    double *u = mxGetPr(prhs[2]);
    double *v = mxGetPr(prhs[3]);
	double *nx = mxGetPr(prhs[4]);
    double *ny = mxGetPr(prhs[5]);
    double *eidM = mxGetPr(prhs[6]);
    double *eidP = mxGetPr(prhs[7]);
    signed char *eidtype = (signed char *)mxGetData(prhs[8]); // int8 ç±»å?

	/* get dimensions */
    size_t Nfp = mxGetM(prhs[7]);
    size_t K = mxGetN(prhs[7]);

	/* allocate output array */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);
	double *dflux  = mxGetPr(plhs[0]);

	/* set number of threads */
    int i,j;
    #ifdef _OPENMP
    #pragma omp parallel for private(j) num_threads(DG_THREADS)
    #endif
	for (i=0;i<K;i++){
        int ind = i*Nfp;
		for(j=0;j<Nfp;j++){
            int iM = (int)eidM[ind]-1; // change index to C type
            int iP = (int)eidP[ind]-1;
            double f_M = h[iM]; // local and adjacent node values
            double varP = h[iP];
            double uM = u[iM], vM = v[iM];
            // double uP = u[iP], vP = v[iP];
            
            // outward normal vector of local element
            double nx_ = nx[ind];
            double ny_ = ny[ind];

            double f_ext; // external values on local nodes
            f_ext = h_ext[iM];

            bc_type type = (bc_type)eidtype[ind];
            // get adjacent values hP, qxP, qyP, considering
            // various boudnary conditions
            double f_P;
            int info = bound_cond(f_M, varP, f_ext, nx_, ny_, type, &f_P);
            // if(info) mexErrMsgTxt("Unknown boundary conditions.");

            double numflux;
			upwind_flux(f_M, f_P, uM, vM, nx_, ny_, &numflux);
            dflux[ind] = numflux;
            ind++;
		}
	}

    return;
}
