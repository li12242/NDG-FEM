#include "mex.h"
#include "mxSWE1d.h"

void evaluateExactAdjacentNodeValue(NdgEdgeType type,
                                    const double nx,
                                    double* fm,
                                    double* fp,
                                    double* fe) {
  switch (type) {
    case NdgEdgeInner:
      break;
    case NdgEdgeSlipWall:
      evaluateSlipWallAdjacentNodeValue(nx, fm, fp);
      break;
    case NdgEdgeNonSlipWall:
      evaluateNonSlipWallAdjacentNodeValue(nx, fm, fp);
      break;
    case NdgEdgeZeroGrad:
      fp[0] = fm[0];
      fp[1] = fm[1];
      break;
    case NdgEdgeClamped:
      fp[0] = 2 * fe[0] - fm[0];
      fp[1] = 2 * fe[1] - fm[1];
      break;
    case NdgEdgeClampedDepth:
      fp[0] = fe[0];
      fp[1] = fm[1];
      break;
    case NdgEdgeClampedVel:
      fp[0] = fm[0];
      fp[1] = 2 * fe[1] - fm[1];
      break;
    case NdgEdgeFlather:
      evaluateFlatherAdjacentNodeValue(nx, fm, fp);
      break;
    default:
      break;
  }
  return;
}

/**
 @brief Implement the hydraostatic reconstruction from Hou et. al. (2013)
 */
void evaluateHydrostaticReconstructValue(double hmin, double* fm, double* fp) {
  double zstar = max(fm[2], fp[2]);
  double um, vm;
  evaluateFlowRateByDeptheThreshold(hmin, fm[0], fm[1], &um);
  evaluateFlowRateByDeptheThreshold(hmin, fp[0], fp[1], &up);

  zstar = min(fm[0] + fm[3], zstar);  // z* = min( \eta^-, z* )
  fm[0] = fm[0] + fm[3] - zstar;
  fp[0] = max(0, fp[0] + fp[3] - zstar) - max(0, fp[3] - zstar);
  fm[1] = fm[0] * um;
  fp[1] = fp[0] * up;
  return;
}

#define NRHS 8
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
  double* eidM = mxGetPr(prhs[2]);
  double* eidP = mxGetPr(prhs[3]);
  double* nx = mxGetPr(prhs[4]);
  signed char* eidtype = (signed char*)mxGetData(prhs[5]);
  double* fphys = mxGetPr(prhs[6]);
  double* fext = mxGetPr(prhs[7]);

  const mwSize* dims = mxGetDimensions(prhs[7]);
  const size_t Np = dims[0];
  const size_t K = dims[1];
  const size_t Nfield = dims[2];

  const size_t TNfp = mxGetM(prhs[2]);
  const size_t ndimOut = 3;
  const mwSize dimOut[3] = {TNfp, K, Nfield};
  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double* he = fext;
  double* hue = he + K * Np;

  double* h = fphys;
  double* hu = h + K * Np;
  double* z = h + 2 * K * Np;

  double* hm = mxGetPr(plhs[0]);
  double* hp = mxGetPr(plhs[1]);
  double* hum = hm + K * TNfp;
  double* hup = hp + K * TNfp;

  for (int k = 0; k < K; k++) {
    for (int n = 0; n < TNfp; n++) {
      const size_t sk = k * TNfp + n;
      NdgEdgeType type = (NdgEdgeType)eidtype[sk];
      if (type == NdgEdgeGaussEdge) {
        continue;
      }
      const int im = (int)eidM[sk] - 1;
      const int ip = (int)eidP[sk] - 1;

      double fm[3] = {h[im], hu[im], z[im]};
      double fp[3] = {h[ip], hu[ip], z[ip]};
      double fe[2] = {he[im], hue[im]};

      const double nx_ = nx[sk];
      evaluateExactAdjacentNodeValue(type, nx_, fm, fp, fe);
      evaluateHydrostaticReconstructValue(hcrit, fm, fp);

      hm[sk] = fm[0];
      hp[sk] = fp[0];
      hum[sk] = fm[1];
      hup[sk] = fp[1];
      hvm[sk] = fm[2];
      hvp[sk] = fp[2];
    }
  }
  return;
}