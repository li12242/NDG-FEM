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
  double* fphys = mxGetPr(prhs[2]);
  double* zgrad = mxGetPr(prhs[3]);

  const mwSize* dims = mxGetDimensions(prhs[2]);
  const size_t Np = dims[0];
  const size_t K = dims[1];

  const size_t NdimOut = 3;
  const mwSize dimOut[3] = {Np, K, NVAR};
  plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double* h = fphys;
  double* z = fphys + 3 * K * Np;
  double* bx = zgrad;
  double* by = zgrad + K * Np;

  double* sourceH = mxGetPr(plhs[0]);
  double* sourceQx = sourceH + K * Np;
  double* sourceQy = sourceQx + K * Np;

  for (int k = 0; k < K; k++) {
    NdgRegionType type = (NdgRegionType)regType[k];
    if (type == NdgRegionWet){

        for (int n = 0; n < Np; n++) {
          int sk = k * Np + n;
          sourceQx[sk] = -gra * (h[sk]) * bx[sk];
          sourceQy[sk] = -gra * (h[sk]) * by[sk];
        }
    }
  }

  return;
}