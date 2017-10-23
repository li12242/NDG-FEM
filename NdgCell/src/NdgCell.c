#include "NdgCell.h"

/** @brief free the memory of NdgCell object. */
void freeNdgCell(NdgCell *cell) {
    free(cell->vr);
    free(cell->vs);
    free(cell->vt);
    free(cell->Nfv);
    freeIntMatrix(cell->FToV);
    free(cell->faceType);
    free(cell->r);
    free(cell->s);
    free(cell->t);
    freeMatrix(cell->M);
    freeMatrix(cell->V);
    freeMatrix(cell->Dr);
    freeMatrix(cell->Ds);
    freeMatrix(cell->Dt);
    free(cell->rq);
    free(cell->sq);
    free(cell->tq);
    freeMatrix(cell->Vq);
    free(cell->Nfp);
    freeIntMatrix(cell->Fmask);
    freeMatrix(cell->LIFT);
    free(cell);
}