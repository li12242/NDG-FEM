#include "nconv2d.h"

void lf_flux(double f_M, double f_P,
             double nx, double ny,
             double *num_flux)
{

    double E_M, G_M, E_P, G_P;
    nodal_flux(f_M, &E_M, &G_M);
    nodal_flux(f_P, &E_P, &G_P);
    double c = max(fabs(f_M), fabs(f_P));

    *num_flux = 0.5 * ((E_M + E_P) * nx + (G_M + G_P) * ny) +
                0.5 * c * (f_M - f_P);
    return;
}

/* @brief calculate the surface flux deviation for strong form.
 * 
 * Usages:
 * 	[dflux] = lf_flux(h, h_ext, nx, ny, eidM, eidP, eidtype);
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    /* check input & output */
    if (nrhs != 7)
        mexErrMsgTxt("Wrong number of input arguments.");

    if (nlhs != 1)
        mexErrMsgTxt("Wrong number of output arguments.");

    /* get inputs */
    double *h = mxGetPr(prhs[0]);
    double *h_ext = mxGetPr(prhs[1]);
    double *nx = mxGetPr(prhs[2]);
    double *ny = mxGetPr(prhs[3]);
    double *eidM = mxGetPr(prhs[4]);
    double *eidP = mxGetPr(prhs[5]);
    signed char *eidtype = (signed char *)mxGetData(prhs[6]); // int8

    /* get dimensions */
    size_t Nfp = mxGetM(prhs[4]);
    size_t K = mxGetN(prhs[4]);
    size_t Np = mxGetM(prhs[0]);

    /* allocate output array */
    plhs[0] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);
    double *dflux = mxGetPr(plhs[0]);

#ifdef _OPENMP /* set number of threads */
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int i = 0; i < K; i++)
    {
        int ind = i * Nfp;
        for (int j = 0; j < Nfp; j++)
        {
            int iM = (int)eidM[ind] - 1; // change index to C type
            int iP = (int)eidP[ind] - 1;
            double f_M = h[iM]; // local and adjacent node values
            double varP = h[iP];

            // outward normal vector of local element
            double nx_ = nx[ind];
            double ny_ = ny[ind];

            double f_ext; // external values on local nodes
            f_ext = h_ext[iM];

            bc_type type = (bc_type)eidtype[ind];
            // get adjacent values, considering various boudnary conditions
            double f_P;
            int info = bound_cond(f_M, varP, f_ext, nx_, ny_, type, &f_P);
            // if(info) mexErrMsgTxt("Unknown boundary conditions.");

            double numflux, E, G;
            lf_flux(f_M, f_P, nx_, ny_, &numflux);
            nodal_flux(f_M, &E, &G);

            dflux[ind] = -numflux + nx_ * E + ny_ * G;
            ind++;
        }
    }

    return;
}