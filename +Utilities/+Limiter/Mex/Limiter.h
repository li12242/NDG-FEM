#ifndef Limiter_H 
#define Limiter_H

#include "mex.h"
#include <math.h>

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

#define real double

/* BasicFunction.c */
void minmod(int n, real *a, real *m);
real sign  (real a);
void FaceMean(int Nfaces, int Nfp, real *h, 
	real *ws, real *sJ, real *Fmask,
	real *face_mean, real *face_len);
void GetLocalVar(int Np, real hmean, real xc, real yc, 
 	real *x, real *y, 
 	real phpx, real phpy, real *h);

/* MatrixUtilities.c */
void MatrixSolver2(real *a, real *f, real *x);

#endif