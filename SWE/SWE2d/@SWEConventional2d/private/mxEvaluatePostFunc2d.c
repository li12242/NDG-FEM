#include "mex.h"
#include "mxSWE2d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 5
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

  /* get inputs */
  double hcrit = mxGetScalar(prhs[0]);
  double* fphys = mxGetPr(prhs[1]);
  double* hc = mxGetPr(prhs[2]);
  double* qxc = mxGetPr(prhs[3]);
  double* qyc = mxGetPr(prhs[4]);

  /* get dimensions */
  const mwSize* dims = mxGetDimensions(prhs[1]);
  const size_t Np = dims[0];
  const size_t K = dims[1];

  const size_t ndimOut = 3;
  const mwSize dimOut[3] = {Np, K, NVAR};
  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double* h = fphys;
  double* qx = fphys + K * Np;
  double* qy = fphys + 2 * K * Np;

  double* h_pos = mxGetPr(plhs[0]);
  double* qx_pos = h_pos + K * Np;
  double* qy_pos = h_pos + 2 * K * Np;

  const double ksi = 0.0;
  // cell area and scalar averages
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
  for (int k = 0; k < K; k++) {
    double hmean = hc[k];
    double qxmean = qxc[k];
    double qymean = qyc[k];

    if (hmean <= ksi) {
      for (int n = 0; n < Np; n++) {
        size_t sk = k * Np + n;
        h[sk] = 0;
        qx[sk] = 0;
        qy[sk] = 0;
      }
      continue;
    }

    double hmin = h[k * Np];
    for (int n = 0; n < Np; n++) {
      hmin = min(hmin, h[k * Np + n]);
    }

    double theta;
    if (hmin < hmean) {
      theta = min((hmean - ksi) / (hmean - hmin), 1.0);
    } else {
      theta = 0.0;
    }

    for (int n = 0; n < Np; n++) {
      size_t sk = k * Np + n;
      h_pos[sk] = theta * (h[sk] - hmean) + hmean;
      qx_pos[sk] = theta * (qx[sk] - qxmean) + qxmean;
      qy_pos[sk] = theta * (qy[sk] - qymean) + qymean;

      if (h_pos[sk] < hcrit) {  // dry nodes
        qx_pos[sk] = 0.0;
        qy_pos[sk] = 0.0;
      }
    }
  }
  return;
}