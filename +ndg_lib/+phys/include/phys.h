#ifndef NDG_FEM_PHYS_H
#define NDG_FEM_PHYS_H

#include "mex.h"
#include "../../+mesh/include/mesh.h"

typedef struct dg_phys
{
    dg_mesh *mesh;
    size_t Nfield;
    double *f_Q;
    double *f_extQ;

} dg_phys;

#endif //NDG_FEM_PHYS_H