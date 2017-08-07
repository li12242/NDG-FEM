#ifndef NDG_FEM_MESH_H
#define NDG_FEM_MESH_H

typedef enum {
    NORMAL = 0,
    SPONGE = 1,
    REFINE = 2,
    WET = 4,
    DRY = 5
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
    int K;
    int Nv;
    int **EToV;
    int *EToR;
    int **EToBS;
    double *vx, *vy, *vz;
    double **J;
    int **EToE, **EToF;
    double **x, **y, **z;
    double **rx, **ry, **rz;
    double **sx, **sy, **sz;
    double **tx, **ty, **tz;
    double *vol;
    double **elen;
    int **eidM;
    int **eidP;
    int **eidtype;
    double **nx, **ny, **nz;
    double **Js;

} dg_mesh;

#endif //NDG_FEM_MESH_H