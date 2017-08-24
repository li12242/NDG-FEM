#ifndef NCONV2D_H
#define NCONV2D_H

#include "mex.h"
#include <math.h>

#define max(a, b) (a > b ? a : b)

typedef enum {
    NORMAL = 0,
    SPONGE = 1,
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

/* nconv2d.c */
int bound_cond(double varM, double varP,
               double f_ext,
               double nx, double ny,
               bc_type type,
               double *f_P);

int nodal_flux(double c, double *E, double *G);

#endif