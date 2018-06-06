#include "conv2d.h"

void flux_term(int Np, int K, double *h,
               double *u, double *v, signed char *EToR,
               double *eflux, double *gflux)
{
    int k, i;
#ifdef _OPENMP
#pragma omp parallel for private(i) num_threads(DG_THREADS)
#endif
    for (k = 0; k < K; k++)
    {
        if ((cell_type)EToR[k] == REFINE)
            continue;

        int ind = k * Np;
        for (i = 0; i < Np; i++)
        {
            nodal_flux(h[ind], u[ind], v[ind],
                       eflux + ind, gflux + ind);
            ind++;
        }
    }
}