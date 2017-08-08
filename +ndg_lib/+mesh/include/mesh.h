#ifndef NDG_FEM_MESH_H
#define NDG_FEM_MESH_H

#include "mex.h"
#include "../../std_cell/include/std_cell.h"

typedef enum {
    NORMAL = 0,
    SPONGE = 1,
    REFINE = 2,
    WET = 4,
    DRY = 5,
    PARWET = 6, // partial wet
} dg_mesh_type;

typedef enum {
    INNER = 0,
    SLIPWALL = 2,
    NONSLIPWALL = 3,
    ZEROGRAD = 4,
    CLAMPED = 5,
    CLAMPEDDEPTH = 6,
    CLAMPEDVEL = 7,
    FLATHER = 8
} dg_bc_type;

typedef struct dg_mesh
{
    dg_cell *cell;
    size_t K;
    size_t Nv;
    double *EToV;
    signed char *EToR;
    signed char *EToBS;
    double *vx, *vy, *vz;
    double *J;
    double *EToE, *EToF;
    double *x, *y, *z;
    double *rx, *ry, *rz;
    double **sx, **sy, *sz;
    double *tx, *ty, *tz;
    double *vol;
    double *elen;
    double *eidM;
    double *eidP;
    signed char *eidtype;
    double *nx, *ny, *nz;
    double *Js;

} dg_mesh;

#endif //NDG_FEM_MESH_H