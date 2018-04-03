#include "mex.h"
#include "mxSWE2d.h"

#define NRHS 4
#define NLHS 1
#define NVAR 3

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  /* check input & output */
  if (nrhs != NRHS) {
    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NRHS);
  }

  if (nlhs != NLHS) {
    mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NLHS);
  }

  double gra = mxGetScalar(prhs[0]);
  signed char* regType = (signed char*)mxGetData(prhs[1]);
  PhysField fphys = convertMexToPhysField(prhs[2]);
  // double* fphys = mxGetPr(prhs[2]);
  double* zgrad = mxGetPr(prhs[3]);

  // const mwSize* dims = mxGetDimensions(prhs[2]);
  // const size_t Np = dims[0];
  // const size_t K = dims[1];
  const size_t Np = fphys.Np;
  const size_t K = fphys.K;

  const size_t NdimOut = 3;
  const mwSize dimOut[3] = {Np, K, NVAR};
  plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  // double* h = fphys;
  // double* z = fphys + 3 * K * Np;
  double* bx = zgrad;
  double* by = zgrad + K * Np;

  PhysField source = convertMexToPhysField(plhs[0]);
  // double* sourceH = mxGetPr(plhs[0]);
  // double* sourceQx = sourceH + K * Np;
  // double* sourceQy = sourceQx + K * Np;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
  for (int k = 0; k < K; k++) {
    NdgRegionType type = (NdgRegionType)regType[k];
    if (type == NdgRegionWet) {
      for (int n = 0; n < Np; n++) {
        int sk = k * Np + n;
        const double h_ = fphys.h[sk];
        source.hu[sk] = -gra * (h_)*bx[sk];
        source.hv[sk] = -gra * (h_)*by[sk];
      }
    }
  }
  return;
}