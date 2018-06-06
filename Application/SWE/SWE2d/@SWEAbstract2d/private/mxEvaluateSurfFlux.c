#include <math.h>
#include "mex.h"
#include "mxSWE2d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

inline void evaluateSurfFluxTerm(const double hmin,  ///< water depth threshold
                                 const double gra,   ///< gravity acceleration
                                 const double h,     ///< water depth
                                 const double hu,    ///< water flux
                                 const double hv,    ///< water flux
                                 double* E,          ///< flux term
                                 double* G           ///< flux term
) {
  double u, v;
  evaluateFlowRateByDeptheThreshold(hmin, h, hu, hv, &u, &v);
  const double huv = h * u * v;
  const double h2 = h * h;
  E[0] = hu;
  G[0] = hv;
  E[1] = h * u * u + 0.5 * gra * h2;
  G[1] = huv;
  E[2] = huv;
  G[2] = h * v * v + 0.5 * gra * h2;
  return;
}

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

  double hcrit = mxGetScalar(prhs[0]);
  double gra = mxGetScalar(prhs[1]);
  double* nx = mxGetPr(prhs[2]);
  double* ny = mxGetPr(prhs[3]);
  // double* fm = mxGetPr(prhs[4]);

  PhysField fm = convertMexToPhysField(prhs[4]);

  const mwSize* dims = mxGetDimensions(prhs[4]);
  const size_t TNfp = fm.Np;
  const size_t K = fm.K;

  const size_t ndimOut = 3;
  const mwSize dimOut[3] = {TNfp, K, NVAR};
  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  PhysField surfFlux = convertMexToPhysField(plhs[0]);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif

  for (int k = 0; k < K; k++) {
    for (int n = 0; n < TNfp; n++) {
      const size_t sk = k * TNfp + n;
      // double fm[4] = {fm.h[sk], hu[sk], hv[sk], z[sk]};
      const double nx_ = nx[sk];
      const double ny_ = ny[sk];

      double E[3], G[3];
      evaluateSurfFluxTerm(hcrit, gra, fm.h[sk], fm.hu[sk], fm.hv[sk], E, G);

      surfFlux.h[sk] = nx_ * E[0] + ny_ * G[0];
      surfFlux.hu[sk] = nx_ * E[1] + ny_ * G[1];
      surfFlux.hv[sk] = nx_ * E[2] + ny_ * G[2];
    }
  }
  return;
}