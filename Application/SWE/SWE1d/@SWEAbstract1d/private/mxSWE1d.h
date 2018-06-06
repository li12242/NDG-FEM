#ifndef MXSWE1D_H
#define MXSWE1D_H

#include "mex.h"
#include <math.h>

#define TAU 1e-6

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

#define EPS 1e-6
#define NVAR 2

typedef enum {
  NdgRegionNormal = 1,
  NdgRegionRefine = 2,
  NdgRegionSponge = 3,
  NdgRegionWet = 4,
  NdgRegionDry = 5,
  NdgRegionPartialWet = 6,
  NdgRegionPartialWetFlood = 7,
  NdgRegionPartialWetDamBreak = 8
} NdgRegionType;

typedef enum {
  NdgEdgeInner = 0,
  NdgEdgeGaussEdge = 1,
  NdgEdgeSlipWall = 2,
  NdgEdgeNonSlipWall = 3,
  NdgEdgeZeroGrad = 4,
  NdgEdgeClamped = 5,
  NdgEdgeClampedDepth = 6,
  NdgEdgeClampedVel = 7,
  NdgEdgeFlather = 8
} NdgEdgeType;

typedef struct {
  size_t Np;     ///< length of 1st dimension
  size_t K;      ///< length of 2nd dimension
  size_t Nfield; ///< length of 3rd dimension
  double *h;
  double *hu;
  double *z;
} PhysField1d;

/** convert mex variable to PhysVolField structure */
inline PhysField1d convertMexToPhysField(const mxArray *mxfield) {
  const mwSize *dims = mxGetDimensions(mxfield);
  PhysField1d field;
  field.Np = dims[0];
  field.K = dims[1];
  field.Nfield = dims[2];
  const size_t Ntmp = field.Np * field.K;

  field.h = mxGetPr(mxfield);
  field.hu = field.h + Ntmp;
  field.z = field.hu + Ntmp;
  return field;
}

/** Evaluate the flow rate depending on the depth threshold */
inline void
evaluateFlowRateByDeptheThreshold(const double hcrit, ///< depth threshold
                                  const double h,     ///< depth
                                  const double hu,    ///< water flux
                                  double *u           ///< result velocity
) {
  if (h > hcrit) {
    //     const double sqrt2 = 1.414213562373095;
    //     double h4 = pow(h, 4);
    //     *u = sqrt2 * h * hu / sqrt( h4 + max( hcrit, h4 ) );
    //     *v = sqrt2 * h * hv / sqrt( h4 + max( hcrit, h4 ) );
    *u = hu / h;
  } else {
    *u = 0.0;
  }

  return;
}

#endif // MXSWE1D_H