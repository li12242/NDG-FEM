#include "mex.h"
#include "mxSWE2d.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef enum {
  NdgEdgeInner = 0,
  NdgEdgeGaussEdge = 1,
  NdgEdgeSlipWall = 2,
  NdgEdgeNonSlipWall = 3,
  NdgEdgeZeroGrad = 4,
  NdgEdgeClamped = 5,
  NdgEdgeClampedDepth = 6,
  NdgEdgeClampedVel = 7,
  NdgEdgeFlather = 8,
  NdgEdgeNonLinearFlather = 9,
  NdgEdgeNonLinearFlatherFlow = 10,
  NdgEdgeNonReflectingFlux = 11
} NdgEdgeType;

void imposeBoundaryCondition(const double gra,  ///< gravity acceleration
                             NdgEdgeType type,  ///< boundary node type
                             const double nx,   ///< outward normal vector
                             const double ny,   ///< outward normal vector
                             const double hM, const double huM,
                             const double hvM, const double zM, const double hE,
                             const double huE, const double hvE,
                             const double zE, double *hP, double *huP,
                             double *hvP, double *zP);

typedef struct {
  int Nfp;     ///< num of face nodes
  int Ne;      ///< num of edges
  int Nfield;  ///< num of field
  double *h;   ///< water depth at local node
  double *hu;  ///< flux at local node
  double *hv;  ///< flux at local node
  double *z;   ///< bottom elevation
} SurfNodeField;

SurfNodeField ConvertMexToSurfField(const mxArray *mxField) {
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

#define NRHS 6
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

  double gra = mxGetScalar(prhs[0]);
  double *nx = mxGetPr(prhs[1]);
  double *ny = mxGetPr(prhs[2]);
  SurfNodeField fM = ConvertMexToSurfField(prhs[3]);
  SurfNodeField fE = ConvertMexToSurfField(prhs[4]);
  signed char *ftype = (signed char *)mxGetData(prhs[5]);

  const size_t ndimOut = 3;
  const mwSize dimOut[3] = {fM.Nfp, fM.Ne, fM.Nfield};
  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
  SurfNodeField fP = ConvertMexToSurfField(plhs[0]);

  const int Ne = fM.Ne;
  const int Nfp = fM.Nfp;

  for (int e = 0; e < Ne; e++) {
    NdgEdgeType type = (NdgEdgeType)ftype[e];  // boundary condition
    for (int n = 0; n < Nfp; n++) {
      const int sk = n + e * Nfp;
      imposeBoundaryCondition(gra, type, nx[sk], ny[sk], fM.h[sk], fM.hu[sk],
                              fM.hv[sk], fM.z[sk], fE.h[sk], fE.hu[sk],
                              fE.hv[sk], fE.z[sk], fP.h + sk, fP.hu + sk,
                              fP.hv + sk, fP.z + sk);
    }
  }
}

/** \brief Impose boundary condition to the next node */
void imposeBoundaryCondition(const double gra,  ///< gravity acceleration
                             NdgEdgeType type,  ///< boundary node type
                             const double nx,   ///< outward normal vector
                             const double ny,   ///< outward normal vector
                             const double hM, const double huM,
                             const double hvM, const double zM, const double hE,
                             const double huE, const double hvE,
                             const double zE, double *hP, double *huP,
                             double *hvP, double *zP) {
  // assign the local node values
  *zP = zM;

  // get next node values
  if (type == NdgEdgeInner) {
    *hP = hM;
    *huP = huM;
    *hvP = hvM;
  } else if (type == NdgEdgeZeroGrad) {
    *hP = hM;
    *huP = huM;
    *hvP = hvM;
  } else if (type == NdgEdgeClamped) {
    *hP = hE;
    *huP = huE;
    *hvP = hvE;
  } else if (type == NdgEdgeClampedDepth) {
    *hP = hE;
    *huP = huM;
    *hvP = hvM;
  } else if (type == NdgEdgeClampedVel) {
    *hP = hM;
    *huP = huE;
    *hvP = hvE;
  } else if (type == NdgEdgeSlipWall) {
    const double qxM = huM;
    const double qyM = hvM;
    double qnM = qxM * nx + qyM * ny;   // outward normal flux
    double qvM = -qxM * ny + qyM * nx;  // outward tangential flux
    // adjacent value
    *hP = hM;
    *huP = (-qnM) * nx - qvM * ny;
    *hvP = (-qnM) * ny + qvM * nx;
  } else if (type == NdgEdgeNonSlipWall) {
    *hP = hM;
    *huP = -huM;
    *hvP = -hvM;
  } else if (type == NdgEdgeFlather) {
    const double uE = huE / hE;
    const double vE = hvE / hE;
    const double unE = uE * nx + vE * ny;   // outward normal flux
    const double uvE = -uE * ny + vE * nx;  // tangential flux
    const double un = unE - sqrt(-gra / *zP) * (hE - hM);
    const double uv = uvE;
    *hP = hM;
    *huP = (un * nx - uv * ny) * hM;
    *hvP = (un * ny + uv * nx) * hM;
  } else if (type == NdgEdgeNonLinearFlather) {
    const double uE = huE / hE;
    const double vE = hvE / hE;
    const double unE = uE * nx + vE * ny;  // outward normal flux
    const double un = huM / hM * nx + hvM / hM * ny;
    const double temp = 0.5 * (un - unE) + sqrt(gra * hE);
    *hP = temp * temp / gra;
    *huP = huE;
    *hvP = hvE;
  } else if (type == NdgEdgeNonLinearFlatherFlow) {
    const double uE = huE / hE;
    const double vE = hvE / hE;
    const double unE = uE * nx + vE * ny;   // outward normal flux
    const double uvE = -uE * ny + vE * nx;  // tangential flux
    const double un = unE - 2 * sqrt(gra * hE) + 2 * sqrt(gra * hM);
    const double uv = uvE;
    *hP = hM;
    *huP = (un * nx - uv * ny) * hM;
    *hvP = (un * ny + uv * nx) * hM;
    //   } else if (type == NdgEdgeNonReflectingFlux) {
    //     const double unM = surf->huM / surf->hM * nx + surf->hvM / surf->hM *
    //     ny; const double RLP = unM + 2 * sqrt(gra * surf->hM); const double
    //     hs = fext->h[idM]; surf->hP = hs; const double RRM = RLP - 4 *
    //     sqrt(gra * hs); const double un = 0.5 * (RRM + RLP); const double uvM
    //     = -surf->huM / surf->hM * ny + surf->hvM / surf->hM * nx; surf->huP =
    //     (un * nx - uvM * ny) * surf->hM; surf->hvP = (un * ny + uvM * nx) *
    //     surf->hM;
  } else {
    mexPrintf("Matlab:%s:Unknown boundary type: %d\n", __FILE__, type);
  }
  return;
}