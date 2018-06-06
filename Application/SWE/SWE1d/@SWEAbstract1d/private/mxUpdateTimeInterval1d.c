#include "mex.h"
#include "mxSWE1d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 6
#define NLHS 1

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

  double hmin = mxGetScalar(prhs[0]);
  double gra = mxGetScalar(prhs[1]);
  double N = mxGetScalar(prhs[2]);
  double *dx = mxGetPr(prhs[3]);
  signed char *regionType = (signed char *)mxGetPr(prhs[4]);
  //   double *fphys = mxGetPr(prhs[5]);

  PhysField1d fphys = convertMexToPhysField(prhs[5]);

  const mwSize *dims = mxGetDimensions(prhs[5]);
  const size_t Np = dims[0];
  const size_t K = dims[1];

  plhs[0] = mxCreateDoubleScalar(0);

  //   double *h = fphys;
  //   double *hu = fphys + K * Np;

  double dt = 1e6;
  for (int k = 0; k < K; k++) {
    NdgRegionType type = (NdgRegionType)regionType[k];
    if (type == NdgRegionDry) {
      continue;
    }
    double dx_ = dx[k];
    for (int n = 0; n < Np; n++) {
      const size_t sk = k * Np + n;
      const double h_ = fphys.h[sk];
      if (h_ < hmin) {
        continue;
      }

      double u;
      evaluateFlowRateByDeptheThreshold(hmin, h_, fphys.hu[sk], &u);
      const double spe = fabs(u);
      const double dtloc = dx_ / (spe + sqrt(gra * h_)) / (2 * N + 1);
      dt = min(dt, dtloc);
    }
  }
  *mxGetPr(plhs[0]) = dt;

  return;
}