#ifndef Limiter_H
#define Limiter_H

#include "mex.h"
#include <math.h>

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

#define real double


void minmod(int n, real *a, real *m);
real sign  (real a);
void faceMean(int Nfaces, int Nfp, int Nfields, real **h,
		real *ws, real *sJ, real *Fmask,
		real *face_mean, real *face_len);
void cellMean(int Np, int Nfields, real **h,
    	real *w, real *J, real *cellMean, real *area);

void getLocalVar(int Np, real hmean, real xc, real yc,
 		real *x, real *y,
 		real phpx, real phpy, real *h);
void meanGradient(int Nsub, real *xv, real *yv, real *hv,
		real xc, real yc, real hc, real *dhdx, real *dhdy);
/* MatrixUtilities.c */
void matrixSolver2(real *a, real *f, real *x);

#endif
