#ifndef STD_CELL_H
#define STD_CELL_H

typedef enum {
    POINT = 0,
    LINE = 1,
    TRI = 2,
    QUAD = 3
} dg_cell_type;

typedef struct dg_cell
{
    dg_cell_type type;
    int Nv;
    double vol;
    double *vr, *vs, *vt;
    int *Nfv;
    int **FToV;
    int Nface;
    dg_cell_type *faceType;
    int N;
    int Np;
    double *r, *s, *t;
    double **V, **M;
    double **Dr, **Ds, **Dt;
    double *w;
    int **Fmask;
    int *Nfp;
    int Nfptotal;
    double **LIFT;
    double *ws;
} dg_cell;

#endif //STD_CELL_H