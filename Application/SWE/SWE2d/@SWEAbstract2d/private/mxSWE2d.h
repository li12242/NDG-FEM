#ifndef __mxSWE_H__
#define __mxSWE_H__

#include "mex.h"
#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

#define EPS 1e-6

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

typedef struct {
  size_t Np;      ///< length of 1st dimension
  size_t K;       ///< length of 2nd dimension
  size_t Nfield;  ///< length of 3rd dimension
  double *h;
  double *hu;
  double *hv;
  double *z;
} PhysField;

/** convert mex variable to PhysVolField structure */
inline PhysField convertMexToPhysField(const mxArray *mxfield) {
  const mwSize *dims = mxGetDimensions(mxfield);
  PhysField field;
  field.Np = dims[0];
  field.K = dims[1];
  field.Nfield = dims[2];
  const size_t Ntmp = field.Np * field.K;

  field.h = mxGetPr(mxfield);
  field.hu = field.h + Ntmp;
  field.hv = field.hu + Ntmp;
  field.z = field.hv + Ntmp;
  return field;
}

/** Evaluate the flow rate depending on the depth threshold */
inline void evaluateFlowRateByDeptheThreshold(
    const double hcrit,  ///< depth threshold
    const double h,      ///< depth
    const double hu,     ///< water flux
    const double hv,     ///< water flux
    double *u,           ///< result velocity
    double *v            ///< velocity
) {
  if (h > hcrit) {
    //     const double sqrt2 = 1.414213562373095;
    //     double h4 = pow(h, 4);
    //     *u = sqrt2 * h * hu / sqrt( h4 + max( hcrit, h4 ) );
    //     *v = sqrt2 * h * hv / sqrt( h4 + max( hcrit, h4 ) );
    *u = hu / h;
    *v = hv / h;
  } else {
    *u = 0.0;
    *v = 0.0;
  }

  return;
}

#endif  //__mxSWE_H__