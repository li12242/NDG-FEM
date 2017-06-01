#ifndef VB_H
#define VB_H

#include "mex.h"
#include <math.h>
#include <omp.h>

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

typedef void (* Wei_Fun_t)(int Nsub, double *grad_x, double *grad_y,
	double *gra_det, double *dfdx, double *dfdy);

void VertLimit(int K, int Np, int Nfaces, int Nfp,
    double *f_v, double *f_max, double *f_min,
	double *fc, double *xc, double *yc,
	double *f_Q, double *x, double *y,
	double *Fmask, double *EToV, double *flim, Wei_Fun_t WeiGrad);

void GetWeiGrad(int Nsub, double *xv, double *yv, double *hv,
    double xc, double yc, double hc,
    double *dhdx, double *dhdy, Wei_Fun_t WeiGrad);

void Grad2Node(int Np, double hmean,
    double xc, double yc,
 	double *x, double *y,
    double phpx, double phpy, double *h);

void MatrixSolver2(double *a, double *f, double *x);

#endif
