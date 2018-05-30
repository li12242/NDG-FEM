#include "mex.h"
#include "mxSWE2d.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
  int Nfp;     ///< num of face nodes
  int Ne;      ///< num of edges
  int Nfield;  ///< num of field
  double *h;   ///< water depth at local node
  double *hu;  ///< flux at local node
  double *hv;  ///< flux at local node
  double *z;   ///< bottom elevation
} SurfNodeField;

//---------------------------------------------------------
/// \brief Implement the hydraostatic reconstruction
/// \details Details about the HR method can be found in Hou et. al. (2013)
void evaluateHydrostaticReconstructValue(
    const double hmin,                     ///< water depth threshold
    const int ind,                         ///< index
    SurfNodeField fM, SurfNodeField fP,    ///< input
    SurfNodeField fM1, SurfNodeField fP1)  ///< result
//---------------------------------------------------------
{
  double zstar = max(fM.z[ind], fP.z[ind]);
  double um, vm, up, vp;
  evaluateFlowRateByDeptheThreshold(hmin, fM.h[ind], fM.hu[ind], fM.hv[ind],
                                    &um, &vm);
  evaluateFlowRateByDeptheThreshold(hmin, fP.h[ind], fP.hu[ind], fP.hv[ind],
                                    &up, &vp);
  const double etaM = fM.h[ind] + fM.z[ind];
  const double etaP = fP.h[ind] + fP.z[ind];
  zstar = min(etaM, zstar);  // z* = min( \eta^-, z* )
  fM1.h[ind] = etaM - zstar;
  fP1.h[ind] = max(0, etaP - zstar) - max(0, fP.z[ind] - zstar);
  fM1.hu[ind] = fM.h[ind] * um;
  fP1.hu[ind] = fP.h[ind] * up;
  fM1.hv[ind] = fM.h[ind] * vm;
  fP1.hv[ind] = fP.h[ind] * vp;
  fM1.z[ind] = zstar;
  fP1.z[ind] = zstar;
  return;
}

//---------------------------------------------------------
SurfNodeField ConvertMexToSurfField(const mxArray *mxField)
//---------------------------------------------------------
{
  SurfNodeField field;
  const mwSize *dims = mxGetDimensions(mxField);  // local phys field dimension
  const int Nfp = dims[0];                        // num of interp nodes
  const int Ne = dims[1];                         // num of elements
  int Nfield;
  if (mxGetNumberOfDimensions(mxField) > 2) {
    Nfield = dims[2];
  } else {
    Nfield = 1;
  }
  field.h = mxGetPr(mxField);
  field.hu = field.h + Nfp * Ne;
  field.hv = field.hu + Nfp * Ne;
  field.z = field.hv + Nfp * Ne;
  field.Nfp = Nfp;
  field.Ne = Ne;
  field.Nfield = Nfield;

  return field;
}

#define NRHS 3
#define NLHS 2

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
  SurfNodeField fM = ConvertMexToSurfField(prhs[1]);
  SurfNodeField fP = ConvertMexToSurfField(prhs[2]);
  const int Nfp = fM.Nfp;
  const int Ne = fM.Ne;

  const size_t ndimOut = 3;
  const mwSize dimOut[3] = {fM.Nfp, fM.Ne, fM.Nfield};
  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
  SurfNodeField fM1 = ConvertMexToSurfField(plhs[0]);
  SurfNodeField fP1 = ConvertMexToSurfField(plhs[1]);

  for (int k = 0; k < Ne; k++) {
    for (int n = 0; n < Nfp; n++) {
      const int sk = k * Nfp + n;
      evaluateHydrostaticReconstructValue(hmin, sk, fM, fP, fM1, fP1);
    }
  }
}