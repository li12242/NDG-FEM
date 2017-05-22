#ifndef __TVB_H
#define __TVB_H

#include "mex.h"
#include <math.h>
#define TOTALERR 1e-10

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )
#define sign(a) ( (a>=0)?1:-1 )

void Grad2Node(int Np, double hmean,
    double xc, double yc,
 	double *x, double *y,
    double phpx, double phpy, double *h);

void MatrixSolver2(double *a, double *f, double *x);
void minmod(int n, double *a, double *m);
double TVB_minmod(double a1, double a2, double dx2, double TVB_factor);

#endif //__TVB_H
