#include "../../@SWEAbstract1d/private/mxSWE1d.h"
#include "mex.h"

#define NRHS 4
#define NLHS 1
#define NVAR 2

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
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
  signed char *regType = (signed char *)mxGetData(prhs[1]);
  PhysField1d fphys = convertMexToPhysField(prhs[2]);
  //   double *fphys = mxGetPr(prhs[2]);
  double *zgrad = mxGetPr(prhs[3]);

  //   const mwSize *dims = mxGetDimensions(prhs[2]);
  const size_t Np = fphys.Np;
  const size_t K = fphys.K;

  const mwSize dimLen[3] = {Np, K, NVAR};
  plhs[0] = mxCreateNumericArray(3, dimLen, mxDOUBLE_CLASS, mxREAL);
  PhysField1d source = convertMexToPhysField(plhs[0]);
  double *bx = zgrad;

  //   double *sourceH = mxGetPr(plhs[0]);
  //   double *sourceQx = sourceH + K * Np;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
  for (int k = 0; k < K; k++) {
    NdgRegionType type = (NdgRegionType)regType[k];
    if (type == NdgRegionDry)
      continue;

    for (int n = 0; n < Np; n++) {
      int sk = k * Np + n;
      source.hu[sk] = -gra * (fphys.h[sk] + fphys.z[sk]) * bx[sk];
    }
  }

  return;
}