#include "../../@SWEAbstract2d/private/mxSWE2d.h"

NdgRegionType getCellWDType(const int Np, const double hmin, double *h,
                            double *z) {
  NdgRegionType type = NdgRegionWet; // initialize to wet element
  int wetNodeNum = 0;
  bool partialWetFlood = true;
  for (int n = 0; n < Np; n++) {
    if (h[n] < hmin) { // dry nodes
      if ((h[n] + z[n] - hmin) > z[n]) {
        partialWetFlood = false;
      }
    } else if (h[n] > hmin) { // wet nodes
      wetNodeNum += 1;
    }
  }

  if (wetNodeNum == 0) { // dry element
    type = NdgRegionDry;
  } else if (wetNodeNum < Np) { // partial wet element
    if (partialWetFlood) {
      type = NdgRegionPartialWetFlood;
    } else {
      type = NdgRegionPartialWetDamBreak;
    }
  }
  return type;
}

#define NRHS 2
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
  PhysField fphys = convertMexToPhysField(prhs[1]);

  const size_t ndimOut = 1;
  const mwSize dimLen[1] = {fphys.K};
  plhs[0] = mxCreateNumericArray(ndimOut, dimLen, mxINT8_CLASS, mxREAL);
  signed char *cellType = (signed char *)mxGetData(plhs[0]);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
  for (int k = 0; k < fphys.K; k++) {
    cellType[k] = (signed char)getCellWDType(
        fphys.Np, hmin, fphys.h + k * fphys.Np, fphys.z + k * fphys.Np);
  }
  return;
}