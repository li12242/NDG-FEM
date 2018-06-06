#include "mex.h"
#include "mxSWE1d.h"

typedef struct {
  double hM;  ///< water depth at local node
  double hP;  ///< water depth at next node
  double huM; ///< flux at local node
  double huP; ///< flux at next node
  double zM;  ///< bottom elevation
  double zP;  ///< bottom elevation
} SurfNodeField;

void imposeBoundaryCondition(NdgEdgeType type, ///< boundary type
                             const double gra, ///< gravity acceleration
                             const double nx,  ///< outward normal vector
                             const int iM,     ///< local boundary node index
                             const int iP, ///< neighbour boundary node index
                             const PhysField1d *fphys, ///< physical field
                             const PhysField1d *fext,  ///< external field
                             SurfNodeField *surf       ///< boundary state
) {
  surf->hM = fphys->h[iM];
  surf->huM = fphys->hu[iM];
  surf->zM = fphys->z[iM];
  surf->zP = fphys->z[iP];
  if (type == NdgEdgeInner) {
    surf->hP = fphys->h[iP];
    surf->huP = fphys->hu[iP];
    surf->zP = fphys->z[iP];
  } else if ((type == NdgEdgeSlipWall) || (type == NdgEdgeNonSlipWall)) {
    surf->hP = surf->hM;
    surf->huP = -surf->huM;
  } else if (type == NdgEdgeZeroGrad) {
    surf->hP = surf->hM;
    surf->huP = surf->huM;
  } else if (type == NdgEdgeClamped) {
    surf->hP = 2 * fext->h[iM] - surf->hM;
    surf->huP = 2 * fext->hu[iM] - surf->huM;
  } else if (type == NdgEdgeClampedDepth) {
    surf->hP = 2 * fext->h[iM] - surf->hM;
    surf->huP = surf->huM;
  } else if (type == NdgEdgeClampedVel) {
    surf->hP = surf->hM;
    surf->huP = 2 * fext->hu[iM] - surf->huM;
  } else if (type == NdgEdgeFlather) {
    const double qn =
        fext->hu[iM] * nx - sqrt(gra * surf->hM) * (fext->h[iM] - surf->hM);

    surf->hP = surf->hM;
    surf->huP = qn * nx;
  }
  return;
}

// void evaluateExactAdjacentNodeValue(NdgEdgeType type, const double nx,
//                                     double *fm, double *fp, double *fe) {
//   switch (type) {
//   case NdgEdgeInner:
//     break;
//   case NdgEdgeSlipWall:
//     evaluateSlipWallAdjacentNodeValue(nx, fm, fp);
//     break;
//   case NdgEdgeNonSlipWall:
//     evaluateNonSlipWallAdjacentNodeValue(nx, fm, fp);
//     break;
//   case NdgEdgeZeroGrad:
//     fp[0] = fm[0];
//     fp[1] = fm[1];
//     break;
//   case NdgEdgeClamped:
//     fp[0] = 2 * fe[0] - fm[0];
//     fp[1] = 2 * fe[1] - fm[1];
//     break;
//   case NdgEdgeClampedDepth:
//     fp[0] = fe[0];
//     fp[1] = fm[1];
//     break;
//   case NdgEdgeClampedVel:
//     fp[0] = fm[0];
//     fp[1] = 2 * fe[1] - fm[1];
//     break;
//   case NdgEdgeFlather:
//     evaluateFlatherAdjacentNodeValue(nx, fm, fp);
//     break;
//   default:
//     break;
//   }
//   return;
// }

// /** \brief Implement the hydraostatic reconstruction from Hou et. al. (2013)
// */ void evaluateHydrostaticReconstructValue(double hmin, double *fm, double
// *fp) {
//   double zstar = max(fm[2], fp[2]);
//   double um, up;
//   evaluateFlowRateByDeptheThreshold(hmin, fm[0], fm[1], &um);
//   evaluateFlowRateByDeptheThreshold(hmin, fp[0], fp[1], &up);

//   zstar = min(fm[0] + fm[2], zstar); // z* = min( \eta^-, z* )
//   fm[0] = fm[0] + fm[2] - zstar;
//   fp[0] = max(0, fp[0] + fp[2] - zstar) - max(0, fp[2] - zstar);
//   fm[1] = fm[0] * um;
//   fp[1] = fp[0] * up;
//   return;
// }

void evaluateHydrostaticReconstructValue(
    const double hmin,  ///< water depth threshold
    SurfNodeField *surf ///< surface node values
) {
  double zstar = max(surf->zM, surf->zP);
  double um, up;
  evaluateFlowRateByDeptheThreshold(hmin, surf->hM, surf->huM, &um);
  evaluateFlowRateByDeptheThreshold(hmin, surf->hP, surf->huP, &up);
  const double etaM = surf->hM + surf->zM;
  const double etaP = surf->hP + surf->zP;
  zstar = min(etaM, zstar); // z* = min( \eta^-, z* )
  surf->hM = etaM - zstar;
  surf->hP = max(0, etaP - zstar) - max(0, surf->zP - zstar);
  surf->huM = surf->hM * um;
  surf->huP = surf->hP * up;
  return;
}

#define NRHS 8
#define NLHS 2
#define NVAR 2

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check input & output */
  if (nrhs != NRHS) {
    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NRHS);
  }

  if (nlhs != NLHS) {
    mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
    mexPrintf("%d output required.\n", NLHS);
  }

  double hcrit = mxGetScalar(prhs[0]);
  double gra = mxGetScalar(prhs[1]);
  double *eidM = mxGetPr(prhs[2]);
  double *eidP = mxGetPr(prhs[3]);
  double *nx = mxGetPr(prhs[4]);
  signed char *eidtype = (signed char *)mxGetData(prhs[5]);
  // double *fphys = mxGetPr(prhs[6]);
  // double *fext = mxGetPr(prhs[7]);
  PhysField1d fphys = convertMexToPhysField(prhs[6]);
  PhysField1d fext = convertMexToPhysField(prhs[7]);

  // const size_t Np = fphys.Np;
  const size_t K = fphys.K;
  // const size_t Nfield = fphys.Nfield;
  // const mwSize *dims = mxGetDimensions(prhs[7]);
  // const size_t Np = dims[0];
  // const size_t K = dims[1];
  // const size_t Nfield = dims[2];

  const size_t TNfp = mxGetM(prhs[2]);
  // const size_t ndimOut = 3;
  const mwSize dimOut[3] = {TNfp, K, NVAR};
  plhs[0] = mxCreateNumericArray(3, dimOut, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(3, dimOut, mxDOUBLE_CLASS, mxREAL);

  PhysField1d fM = convertMexToPhysField(plhs[0]);
  PhysField1d fP = convertMexToPhysField(plhs[1]);
  // double *he = fext;
  // double *hue = he + K * Np;

  // double *h = fphys;
  // double *hu = h + K * Np;
  // double *z = h + 2 * K * Np;

  // double *hm = mxGetPr(plhs[0]);
  // double *hp = mxGetPr(plhs[1]);
  // double *hum = hm + K * TNfp;
  // double *hup = hp + K * TNfp;

  for (int k = 0; k < K; k++) {
    for (int n = 0; n < TNfp; n++) {
      const size_t sk = k * TNfp + n;
      NdgEdgeType type = (NdgEdgeType)eidtype[sk];
      if (type == NdgEdgeGaussEdge) {
        continue;
      }
      const int im = (int)eidM[sk] - 1;
      const int ip = (int)eidP[sk] - 1;

      // double fm[3] = {h[im], hu[im], z[im]};
      // double fp[3] = {h[ip], hu[ip], z[ip]};
      // double fe[2] = {he[im], hue[im]};
      SurfNodeField surf;
      const double nx_ = nx[sk];
      imposeBoundaryCondition(type, gra, nx_, im, ip, &fphys, &fext, &surf);
      // evaluateExactAdjacentNodeValue(type, nx_, fm, fp, fe);
      evaluateHydrostaticReconstructValue(hcrit, &surf);
      // evaluateHydrostaticReconstructValue(hcrit, fm, fp);

      // hm[sk] = fm[0];
      // hp[sk] = fp[0];
      // hum[sk] = fm[1];
      // hup[sk] = fp[1];
      fM.h[sk] = surf.hM;
      fM.hu[sk] = surf.huM;
      fP.h[sk] = surf.hP;
      fP.hu[sk] = surf.huP;
    }
  }
  return;
}