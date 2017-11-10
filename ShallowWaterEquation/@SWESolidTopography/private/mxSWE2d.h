#ifndef __mxSWE_H__
#define __mxSWE_H__

#define Nfield 6
#define Nvar 3

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

#define EPS 1e-6

typedef enum {
  NdgRegionNormal = 1,
  NdgRegionRefine = 2,
  NdgRegionSponge = 3,
  NdgRegionWet = 4,
  NdgRegionDry = 5,
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
  NdgEdgeFlather = 8,
} NdgEdgeType;

/**
 Evaluate the flow rate depending on the depth threshold
 */
void evaluateFlowRateByDeptheThreshold(const double hcrit,
                                       const double h, const double hu, const double hv,
                                       double *u, double *v);
/**
 Evaluate the flow rate depending on the cell states
 */
void evaluateFlowRateByCellState(const NdgRegionType type,
                                 const double h, const double hu,
                                 const double hv, double *u, double *v);

void evaluateSlipWallAdjacentNodeValue(const double nx, const double ny, double *fm, double *fp);

void evaluateNonSlipWallAdjacentNodeValue(const double nx, const double ny, double *fm, double *fp);

void evaluateFlatherAdjacentNodeValue(double nx, double ny, double *fm, double *fe);

#endif //__mxSWE_H__