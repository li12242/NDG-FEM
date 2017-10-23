#include "NdgCell.h"
#include "blas.h"

int testDerivativeMatrix(NdgCell *cell);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
    NdgCell *cell = mxGetNdgCell(prhs[0]);

    mxPrintNdgCellInfo(cell, "triangle");
    testDerivativeMatrix(cell);
    freeNdgCell(cell);
    return;
}

int testDerivativeMatrix(NdgCell *cell){
    
    ptrdiff_t Np = (size_t)cell->Np;
    ptrdiff_t oneInt = 1;
    char *chn = "N";
    /* scalar values to use in dgemm */
    double one = 1.0, zero = 0.0;
    
    double *C = (double *)calloc(cell->Np, sizeof(double));
    double *dr = cell->Dr[0];
    double *ds = cell->Ds[0];
    double *dt = cell->Dt[0];
    double *r = cell->r;
    double *s = cell->s;
    double *t = cell->t;
    
    dgemm(chn, chn, &Np, &oneInt, &Np, &one, dr, &Np, r, &Np, &zero, C, &Np);

    mexPrintf("Dr * r = \n");
    for (int n=0; n<Np; n++) {
        mexPrintf("%f ", C[n]);
    }
    mexPrintf("\n");
    
    dgemm(chn, chn, &Np, &oneInt, &Np, &one, dr, &Np, s, &Np, &zero, C, &Np);
    
    mexPrintf("Dr * s = \n");
    for (int n=0; n<Np; n++) {
        mexPrintf("%f ", C[n]);
    }
    mexPrintf("\n");
    return 0;
}