#include "NdgMesh.h"



/** Free the memory of NdgMeshUnion structure. */
void
freeNdgUMeshUnion(NdgUMeshUnion *mesh){

    freeNdgCell(mesh->cell);
    freeIntMatrix(mesh->EToV);
    freeIntMatrix(mesh->EToE);
    freeIntMatrix(mesh->EToF);
    freeIntMatrix(mesh->EToM);
    free(mesh->EToB[0]);
    free(mesh->EToB);
    free(mesh->EToR);
    free(mesh->vx);
    free(mesh->vy);
    free(mesh->vz);
    freeMatrix(mesh->x);
    freeMatrix(mesh->y);
    freeMatrix(mesh->z);
    freeMatrix(mesh->J);
    freeMatrix(mesh->rx);
    freeMatrix(mesh->ry);
    freeMatrix(mesh->rz);
    freeMatrix(mesh->sx);
    freeMatrix(mesh->sy);
    freeMatrix(mesh->sz);
    freeMatrix(mesh->tx);
    freeMatrix(mesh->ty);
    freeMatrix(mesh->tz);
    freeIntMatrix(mesh->eidM);
    freeIntMatrix(mesh->eidP);
    free(mesh->eidType[0]);
    free(mesh->eidType);
    freeMatrix(mesh->nx);
    freeMatrix(mesh->ny);
    freeMatrix(mesh->nz);
    freeMatrix(mesh->Js);
    
    free(mesh);
}

