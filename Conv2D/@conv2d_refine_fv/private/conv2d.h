#ifndef CONV2D_H
#define CONV2D_H

#include "mex.h"
// #include "../../../+ndg_lib/+phys/include/phys.h"

typedef enum {
    NORMAL = 0,
    SPONGE = 1,
    REFINE = 2,
    WET = 4,
    DRY = 5,
    PARWET = 6, // partial wet
} cell_type;

typedef enum {
    Inner = 0,
    SlipWall = 2,
    NSlipWall = 3,
    ZeroGrad = 4,
    Clamped = 5,
    ClampedDepth = 6,
    ClampedVel = 7,
    Flather = 8,
} bc_type;

/* conv2d.c */
int nodal_flux(double c, double u, double v, double *E, double *G);
int bound_cond(double varM, double varP, double f_ext,
               double nx, double ny, bc_type type, double *f_P);

void dvecm(double N, double alpha, double *x, double *y, double *t);
void dvecd(double N, double alpha, double *x, double *y, double *t);
void upwind_flux(double f_M, double f_P, double uM, double vM,
                 double nx, double ny, double *numflux);
#endif //CONV2D_H