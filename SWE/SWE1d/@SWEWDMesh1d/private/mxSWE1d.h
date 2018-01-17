#ifndef MXSWE1D_H
#define MXSWE1D_H

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
  NdgRegionPartialWD = 6,
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
  NdgEdgeRefine = 9,
} NdgEdgeType;

void evaluateFlowRateByDeptheThreshold(double hcrit,
                                       double h,
                                       double hu,
                                       double* u);

void evaluateFlowRateByCellState(const NdgRegionType type,
                                 const double h,
                                 const double hu,
                                 double* u);

void evaluateSurfFluxTerm(double hmin,
                          double gra,
                          double h,
                          double hu,
                          double* E);

void evaluateSlipWallAdjacentNodeValue(const double nx, double* fm, double* fp);

void evaluateNonSlipWallAdjacentNodeValue(const double nx,
                                          double* fm,
                                          double* fp);

void evaluateFlatherAdjacentNodeValue(double nx, double* fm, double* fe);

#endif  // MXSWE1D_H