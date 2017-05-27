#ifndef _SWE_H
#define _SWE_H

#include "mex.h"
#include <math.h>

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

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

void nodal_flux(double hcrit, double gra,
    double h, double qx, double qy, double z,
    double *Eh, double *Eqx, double *Eqy,
    double *Gh, double *Gqx, double *Gqy);

int bound_cond(double *varM, double *varP, double *f_ext,
    double nx, double ny, bc_type type, double *f_P);

#endif
