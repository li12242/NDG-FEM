#include "../../@SWEAbstract1d/private/mxSWE1d.h"
#include "mex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void evaluateNodalFlux(double hcrit, double gra, double h, double qx, double z,
                       double *Eh, double *Eqx) {

  double h2 = 0.5 * gra * (h * h - z * z);
  if (h > hcrit) {

    *Eh = qx;
    *Eqx = (qx * qx / h + h2);
  } else { // for dry nodes
    *Eh = 0;
    *Eqx = 0;
  }
  return;
}

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

  double hcrit = mxGetScalar(prhs[0]);
  double gra = mxGetScalar(prhs[1]);
  signed char *regType = (signed char *)mxGetData(prhs[2]);
  //   double *fphys = mxGetPr(prhs[3]);
  PhysField1d fphys = convertMexToPhysField(prhs[3]);
  const size_t Np = fphys.Np;
  const size_t K = fphys.K;
  //   const mwSize *dims = mxGetDimensions(prhs[3]);

  //   const size_t Np = dims[0];
  //   const size_t K = dims[1];

  //   const size_t ndimOut = 3;
  const mwSize dimOut[3] = {Np, K, NVAR};
  plhs[0] = mxCreateNumericArray(3, dimOut, mxDOUBLE_CLASS, mxREAL);
  PhysField1d E = convertMexToPhysField(plhs[0]);

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
      evaluateNodalFlux(hcrit, gra, fphys.h[sk], fphys.hu[sk], fphys.z[sk],
                        E.h + sk, E.hu + sk);
    }
  }

  return;
}