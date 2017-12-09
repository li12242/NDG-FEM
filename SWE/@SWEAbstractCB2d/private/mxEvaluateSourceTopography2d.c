#include "mex.h"
#include "mxSWE2d.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double gra = mxGetScalar(prhs[0]);
    signed char *regType = (signed char *)mxGetData(prhs[1]);
    double *fphys = mxGetPr(prhs[2]);
    
    const mwSize *dims = mxGetDimensions(prhs[2]);
    
    const size_t Np = dims[0];
    const size_t K = dims[1];
    
    const size_t ndimOut = 3;
    const mwSize dimOut[3] = {Np, K, Nvar};
    plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
    
    double *h = fphys;
    double *bot = fphys + 3*K*Np;
    double *bx = bot + K*Np;
    double *by = bot + 2*K*Np;
    
    
    double *sourceH = mxGetPr(plhs[0]);
    double *sourceQx = sourceH + K*Np;
    double *sourceQy = sourceQx + K*Np;
    
    for (int k = 0; k<K; k++) {
        NdgRegionType type = (NdgRegionType) regType[k];
        if (type == NdgRegionDry)
            continue;

        for (int n = 0; n<Np; n++) {
            int sk = k*Np + n;
            sourceQx[sk] = -gra * ( h[sk] + bot[sk] ) * bx[sk];
            sourceQy[sk] = -gra * ( h[sk] + bot[sk] ) * by[sk];
        }
    }
    
    return;
}