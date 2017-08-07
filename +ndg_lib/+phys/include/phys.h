#ifndef NDG_FEM_PHYS_H
#define NDG_FEM_PHYS_H

typedef struct dg_phys
{
    int Nfield;
    double *f_Q;
    double *f_extQ;

} dg_phys;

#endif //NDG_FEM_PHYS_H