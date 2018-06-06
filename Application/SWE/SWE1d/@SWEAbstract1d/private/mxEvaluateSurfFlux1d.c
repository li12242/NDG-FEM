#include "mex.h"
#include "mxSWE1d.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void evaluateSurfFluxTerm(double hmin, ///< wet-dry depth threshold
                          double gra,  ///< gravity acceleration
                          double h,    ///< water depth
                          double hu,   ///< flux
                          double *Eh,  ///< flux term
                          double *Ehu  ///< flux term
) {
  double u;
  evaluateFlowRateByDeptheThreshold(hmin, h, hu, &u);
  const double h2 = h * h;
  *Eh = hu;
  *Ehu = h * u * u + 0.5 * gra * h2;
  return;
}

#define NRHS 4
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

  double hcrit = mxGetScalar(prhs[0]);
  double gra = mxGetScalar(prhs[1]);
  double *nx = mxGetPr(prhs[2]);
  // double* fm = mxGetPr(prhs[3]);
  PhysField1d fm = convertMexToPhysField(prhs[3]);
  const size_t TNfp = fm.Np;
  const size_t K = fm.K;
  // const mwSize *dims = mxGetDimensions(prhs[3]);
  // const size_t TNfp = dims[0];
  // const size_t K = dims[1];

  // const size_t ndimOut = 3;
  const mwSize dimOut[3] = {TNfp, K, NVAR};
  plhs[0] = mxCreateNumericArray(3, dimOut, mxDOUBLE_CLASS, mxREAL);

  PhysField1d flux = convertMexToPhysField(plhs[0]);
  // double *Fhs = mxGetPr(plhs[0]);
  // double *Fqxs = Fhs + TNfp * K;

  // double *h = fm;
  // double *hu = fm + K * TNfp;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif

  for (int k = 0; k < K; k++) {
    for (int n = 0; n < TNfp; n++) {
      const size_t sk = k * TNfp + n;
      const double nx_ = nx[sk];

      evaluateSurfFluxTerm(hcrit, gra, fm.h[sk], fm.hu[sk], flux.h + sk,
                           flux.hu + sk);

      flux.h[sk] *= nx_;
      flux.hu[sk] *= nx_;
    }
  }
  return;
}