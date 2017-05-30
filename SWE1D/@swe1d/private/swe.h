#ifndef SWE_H
#define SWE_H

#include "mex.h"
#include <math.h>

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

typedef enum {
    NORMAL = 0,
    SPONGE = 1,
    WET = 4,
    DRY = 5,
} cell_type;

typedef enum {
    INNER = 0,
    SLIPWALL = 2,
    NSLIPWALL = 3,
    ZEROGRAD = 4,
    CLAMP = 5,
} bc_type;

void nodal_flux(double hcrit, double gra,
    double h, double qx, double z, double *Fh, double *Fq);

void reduce_nodal_flux(double hcrit, double gra,
    double h, double qx, double *Fh, double *Fq);
#endif
