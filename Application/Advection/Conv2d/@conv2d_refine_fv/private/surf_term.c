#include "mex.h"
#include "conv2d.h"

#define DEBUG 0

void surf_term(int Nfp, int K, double *h, double *h_ext,
               double *u, double *v, double *nx, double *ny,
               double *eidM, double *eidP, signed char *eidtype,
               signed char *EToR,
               double *dflux)
{
    int i, j;
#ifdef _OPENMP
#pragma omp parallel for private(j) num_threads(DG_THREADS)
#endif
    for (i = 0; i < K; i++)
    {
        if ((cell_type)EToR[i] == REFINE)
            continue;
        int ind = i * Nfp;
        for (j = 0; j < Nfp; j++)
        {
            int iM = (int)eidM[ind] - 1; // change index to C type
            int iP = (int)eidP[ind] - 1;
            double f_M = h[iM]; // local and adjacent node values
            double hP = h[iP];
            double uM = u[iM], vM = v[iM];

            // outward normal vector of local element
            double nx_ = nx[ind];
            double ny_ = ny[ind];

            double f_ext; // external values on local nodes
            f_ext = h_ext[iM];

            bc_type type = (bc_type)eidtype[ind];
            // get adjacent values hP, qxP, qyP, considering
            // various boudnary conditions
            double f_P;
            int info = bound_cond(f_M, hP, f_ext, nx_, ny_, type, &f_P);
            // if(info) mexErrMsgTxt("Unknown boundary conditions.");

            double numflux, E, G;
            upwind_flux(f_M, f_P, uM, vM, nx_, ny_, &numflux);
            nodal_flux(f_M, uM, vM, &E, &G);

#if DEBUG
            mexPrintf("n = %d, k = %d, num_flux = %e, E = %f, G = %f\n",
                      j, i, numflux, E, G);
#endif
            dflux[ind] = -numflux + nx_ * E + ny_ * G;
            ind++;
        }
    }
    return;
}
