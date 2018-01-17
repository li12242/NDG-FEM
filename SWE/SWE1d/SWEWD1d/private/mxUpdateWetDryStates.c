#include "mex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 3
#define NLHS 1

typedef enum {
    NdgRegionNormal = 1,
    NdgRegionRefine = 2,
    NdgRegionSponge = 3,
    NdgRegionWet = 4,
    NdgRegionDry = 5,
    NdgRegionPartialWD = 6,
} NdgRegionType;

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

/**
 *
 * [ EToR ] = mxUpdateWetDryStates( hmin, EToR, fphys )
 */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /* check input & output */
    if (nrhs != NRHS) {
        mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
        mexPrintf("%d inputs required.\n", NRHS);
    }

    if (nlhs != NLHS) {
        mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
        mexPrintf("%d inputs required.\n", NLHS);
    }

    double hmin = mxGetScalar(prhs[0]);
    signed char* regionType = (signed char*)mxGetPr(prhs[1]);
    double* fphys = mxGetPr(prhs[2]);

    const mwSize* dims = mxGetDimensions(prhs[2]);
    const size_t Np = dims[0];
    const size_t K = dims[1];

    plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
    double* theta = mxGetPr(plhs[0]);
    double* h = fphys;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k = 0; k < K; k++) {
        for (int n = 0; n < Np; n++) {
            const int sk = k * Np + n;
            theta[sk] = 0;
        }

        for (int n = 0; n < Np; n++) {
            const int sk = k * Np + n;

            if (h[sk] >= hmin) {
                const int skm1 = k * Np + max(0, n - 1);
                const int skp1 = k * Np + min(Np - 1, n + 1);
                theta[sk] = 1;
                theta[skm1] = 1;
                theta[skp1] = 1;
            }
        }
    }

    return;
}