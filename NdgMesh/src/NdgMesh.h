#ifndef UMESH_H
#define UMESH_H

#include "NdgCell.h"
#include "NdgNetCDF.h"

/** enumeration for the physical cell types */
typedef enum {
    NdgRegionNromal = 1,
    NdgRegionRefine = 2,
    NdgRegionSponge = 3,
    NdgRegionWet = 4, // SWE solver
    NdgRegionDry = 5, // SWE solver
} NdgRegionType;

/** enumeration for the edge types */
typedef enum {
    NdgEdgeInner = 0, ///< inner edge
    NdgEdgeGauss = 1, ///< gauss integral edges
    NdgEdgeSlipWall = 2, ///< slip wall edges
    NdgEdgeNSlipWall = 3, ///< non-slip wall edges
    NdgEdgeZeroGrad = 4, ///< open boundary condition for zero gradient
    NdgEdgeClamped = 5,///< open boundary condition for clamped condition
} NdgEdgeType;

/** unstructed mesh object of single type standard cell */
typedef struct NdgUMeshUnion {
    
    NdgCell *cell; ///< standard cell object
    int K; ///< total number of local cell
    int Nv; ///< total number of vertex
    int **EToV; ///< vertex index in each cell
    int **EToE; ///< adjacent cell index
    int **EToF; ///< face index of adjacent cell
    int **EToM; ///< mesh index of adjacent cell
    NdgEdgeType **EToB; ///< edge type of each face
    NdgRegionType *EToR; ///< region type of each cell
    double *vx, *vy, *vz; ///< verex coordinate
    double **x, **y, **z; ///< coordinate of interpolation nodes
    double **J; ///< determination of Jacobian matrix
    double **rx, **ry, **rz;
    double **sx, **sy, **sz;
    double **tx, **ty, **tz;
    
    int **eidM, **eidP;
    NdgEdgeType **eidType;
    double **nx, **ny, **nz; ///< outward normal vector
    double **Js;
    
} NdgUMeshUnion;

/** free the memory of NdgMeshUnion structure. */
void freeNdgUMeshUnion(NdgUMeshUnion *mesh);

#include "NdgMeshIO.h"

#endif