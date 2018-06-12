#include <math.h>
#include "mex.h"
#include "mxSWE2d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/** */
void evaluateNodalFlux(double hcrit,  ///< water depth threshold
                       double gra,    ///< gravity accerelation
                       double h,      ///< water depth
                       double qx,     ///< flux
                       double qy,     ///< flux
                       double z,      ///< bottom elevation
                       double *Eh,    ///< flux term
                       double *Eqx,   ///< flux term
                       double *Eqy,   ///< flux term
                       double *Gh,    ///< flux term
                       double *Gqx,   ///< flux term
                       double *Gqy    ///< flux term
) {
  if (h > hcrit) {
    double h2 = 0.5 * gra * (h * h - z * z);
    double huv = qx * qy / h;
    *Eh = qx;
    *Gh = qy;
    *Eqx = (qx * qx / h + h2);
    *Gqx = huv;
    *Eqy = huv;
    *Gqy = (qy * qy / h + h2);
  } else {  // for dry nodes
    *Eh = 0;
    *Eqx = 0;
    *Eqy = 0;
    *Gh = 0;
    *Gqx = 0;
    *Gqy = 0;
  }
  return;
}

#define NRHS 4
#define NLHS 2
#define NVAR 3

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
  signed char *regType = (signed char *)mxGetData(prhs[2]);
  double *fphys = mxGetPr(prhs[3]);

  const mwSize *dims = mxGetDimensions(prhs[3]);

  const size_t Np = dims[0];
  const size_t K = dims[1];

  const size_t ndimOut = 3;
  const mwSize dimOut[3] = {Np, K, NVAR};
  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double *E = mxGetPr(plhs[0]);
  double *G = mxGetPr(plhs[1]);

  double *Eh = E;
  double *Ehu = E + K * Np;
  double *Ehv = E + 2 * K * Np;
  double *Gh = G;
  double *Ghu = G + K * Np;
  double *Ghv = G + 2 * K * Np;

  double *h = fphys;
  double *hu = fphys + K * Np;
  double *hv = fphys + 2 * K * Np;
  double *z = fphys + 3 * K * Np;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
  for (int k = 0; k < K; k++) {
    NdgRegionType type = (NdgRegionType)regType[k];
    if (type == NdgRegionDry) {
      continue;
    }

    for (int n = 0; n < Np; n++) {
      size_t sk = k * Np + n;
      evaluateNodalFlux(hcrit, gra, h[sk], hu[sk], hv[sk], z[sk], Eh + sk,
                        Ehu + sk, Ehv + sk, Gh + sk, Ghu + sk, Ghv + sk);
    }
  }

  return;
}