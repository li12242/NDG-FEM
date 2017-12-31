#include <math.h>
#include "mex.h"
#include "mxSWE2d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define NRHS 4
#define NLHS 1

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

  double hcrit = mxGetScalar(prhs[0]);
  double gra = mxGetScalar(prhs[1]);
  double* nx = mxGetPr(prhs[2]);
  double* fm = mxGetPr(prhs[3]);

  const mwSize* dims = mxGetDimensions(prhs[4]);
  const size_t TNfp = dims[0];
  const size_t K = dims[1];

  const size_t ndimOut = 3;
  const mwSize dimOut[3] = {TNfp, K, NVAR};
  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double* Fhs = mxGetPr(plhs[0]);
  double* Fqxs = Fhs + TNfp * K;
  double* Fqys = Fhs + 2 * TNfp * K;

  double* h = fm;
  double* hu = fm + K * TNfp;
  double* hv = fm + 2 * K * TNfp;
  double* z = fm + 3 * K * TNfp;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif

  for (int k = 0; k < K; k++) {
    for (int n = 0; n < TNfp; n++) {
      const size_t sk = k * TNfp + n;
      double fm[4] = {h[sk], hu[sk], hv[sk], z[sk]};
      const double nx_ = nx[sk];
      const double ny_ = ny[sk];

      double E[3], G[3];
      evaluateSurfFluxTerm(hcrit, gra, fm[0], fm[1], fm[2], E, G);

      Fhs[sk] = nx_ * E[0] + ny_ * G[0];
      Fqxs[sk] = nx_ * E[1] + ny_ * G[1];
      Fqys[sk] = nx_ * E[2] + ny_ * G[2];
    }
  }
  return;
}