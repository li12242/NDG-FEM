#include "mex.h"
#include "mxSWE2d.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* check input & output */
    if (nrhs != 5){
        mexErrMsgIdAndTxt("Matlab:mxUpdateTimeInterval2d:InvalidNumberInput",
                          "5 inputs required.");
    }
    
    if (nlhs != 1){
        mexErrMsgIdAndTxt("Matlab:mxUpdateTimeInterval2d:InvalidNumberOutput",
                          "2 output required.");
    }
    
    double gra = mxGetScalar(prhs[0]);
    double N = mxGetScalar(prhs[1]);
    double *dx = mxGetPr(prhs[2]);
    signed char *regionType = (signed char *)mxGetPr(prhs[3]);
    double *fphys = mxGetPr(prhs[4]);
    
    const mwSize *dims = mxGetDimensions(prhs[4]);
    const size_t Np = dims[0];
    const size_t K = dims[1];
    
    plhs[0] = mxCreateDoubleScalar(0);
    
    double *h = fphys;
    double *hu = fphys + K*Np;
    double *hv = fphys + 2*K*Np;
    
    double dt = 1e6;
    for (int k=0; k<K; k++) {
        NdgRegionType type = (NdgRegionType)regionType[k];
        if ( type == NdgRegionDry ) {
            continue;
        }
        double dx_ = dx[k];
        for (int n=0; n<Np; n++) {
            const size_t sk = k*Np + n;
            double u, v;
            evaluateFlowRateByCellState(type, h[sk], hu[sk], hv[sk], &u, &v);
            const double spe = sqrt( u*u + v*v );
            const double dtloc =  dx_ /( spe + sqrt(gra * h[sk]) )/( 2*N+1 );
            dt = min(dt, dtloc);
        }
    }
    *mxGetPr( plhs[0] ) = dt;
    
    return;
}