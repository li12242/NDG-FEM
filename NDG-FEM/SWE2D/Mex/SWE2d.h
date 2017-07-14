#ifndef SWE2D_H
#define SWE2D_H

// header
#include "mex.h"
#include <math.h>

// real type
#define real double

// operations
#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

// subroutine
void SWE_NodalFlux2d(real hcrit, real gra,
                     real h, real qx, real qy,
                     real *Eh, real *Eqx, real *Eqy,
                     real *Gh, real *Gqx, real *Gqy);

void SWE_HLL2d(real hmin, real gra, real hM, real hP,
               real qnM, real qnP, real qvM, real qvP,
               real *Fhn, real *Fqxn, real *Fqyn);

#endif
