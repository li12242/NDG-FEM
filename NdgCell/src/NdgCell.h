#ifndef STD_CELL_H
#define STD_CELL_H

#include "Utils.h"

/** Enumerations for standard cell type */
typedef enum {
    NdgPoint = 0, ///< Enumerations for point cell
    NdgLine = 1, ///< Enumerations for line cell
    NdgTri = 2, ///< Enumerations for triangular cell
    NdgQuad = 3, ///< Enumerations for quadrilateral cell
    NdgPrismTri = 4, ///< Enumerations for triangular prism cell
    NdgPrismQuad = 5, ///< Enumerations for quadrilateral prism cell
} NdgCellType;

/** Standard cell structure */
typedef struct NdgCell {
    
    /** Maximum order of basis function */
    int N;
    /** Enumeration for cell type */
    NdgCellType type;
    /** Number of vertice */
    int Nv;
    /** Vertice coordinate */
    double *vr, *vs, *vt;
    /**  */
    double vol;
    /** Number of vertice on each face */
    int *Nfv;
    /** Vertice index on each face */
    int **FToV;
    /** Number of faces */
    int Nface;
    /**  */
    NdgCellType *faceType;
    /** Number of interpolation points */
    int Np;
    /** Coordinates for interpolation points */
    double *r, *s, *t;
    int **Fmask;
    int *Nfp;
    int TNfp;
    double **V;
    double **M;
    double **Dr; ///< Matrix of derivative basis function values at interpolation nodes
    double **Ds; ///< Matrix of derivative basis function values at interpolation nodes
    double **Dt;
    double **Drt; ///< transpose of Dr
    double **Dst; ///< matrix Ds and its transpose
    double **Dtt; ///< matrix Dt and its transpose
    double **LIFT; ///<
    double **LIFTt; ///<
    /**  */
    double **Vq;
    /** Number of quadrature points */
    int Nq;
    /** Coordinates for quadrature points */
    double *rq, *sq, *tq;
    /** Integral weights for each quadrature points */
    double *wq;
    
} NdgCell;

/** free the memory of NdgCell structure. */
void freeNdgCell(NdgCell *cell);

#include "NdgCellIO.h"

#endif //STD_CELL_H