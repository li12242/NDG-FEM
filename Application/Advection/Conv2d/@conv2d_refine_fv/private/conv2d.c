#include "conv2d.h"

int nodal_flux(double c, double u, double v,
               double *E, double *G)
{
    *E = u * c;
    *G = v * c;
    return 0;
}

void upwind_flux(double f_M, double f_P, double uM, double vM,
                 double nx, double ny, double *numflux)
{

    const double unM = uM * nx + vM * ny;
    if (unM > 0)
    {
        *numflux = f_M * unM;
    }
    else
    {
        *numflux = f_P * unM;
    }
    return;
}

int bound_cond(double varM, double varP, double f_ext,
               double nx, double ny, bc_type type,
               double *f_P)
{

    switch (type)
    {
    case Inner:
    case SlipWall:
    case NSlipWall:
        *f_P = varP;
        break;
    case ZeroGrad:
        *f_P = varM;
        break;
    case Clamped:
        *f_P = f_ext;
        break;
    case ClampedDepth:
    case ClampedVel:
    case Flather:
        *f_P = 2 * f_ext - varM;
        break;
    }
    return 0;
}

/* 
 * @brief double precision vector multiply operation. 
 * t = x .* y 
 */
void dvecm(double N, double alpha, double *x, double *y, double *t)
{
    int i;
    // #ifdef _OPENMP
    // #pragma omp parallel for num_threads(DG_THREADS)
    // #endif
    for (i = 0; i < N; i++)
    {
        t[i] = t[i] + alpha * x[i] * y[i];
    }
}

/* 
 * @brief double precision vector divide operation. 
 * t = x .* y 
 */
void dvecd(double N, double alpha, double *x, double *y, double *t)
{
    int i;
    // #ifdef _OPENMP
    // #pragma omp parallel for num_threads(DG_THREADS)
    // #endif
    for (i = 0; i < N; i++)
    {
        t[i] = t[i] + alpha * x[i] / y[i];
    }
}