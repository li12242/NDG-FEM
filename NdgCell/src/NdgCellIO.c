//
//  NdgCellIO.c
//  NDG-FEM
//
//  Created by li12242 on 17/10/21.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#include "NdgCellIO.h"

/** @brief allocate memory for NdgCell object. */
static NdgCell *mallocNdgCell() {
    NdgCell *cell = (NdgCell *)calloc(1, sizeof(NdgCell));
    return cell;
}

/** @brief get NdgCellType vector from Matlab class variable. */
static NdgCellType *mxGetPropertyFaceTypeVector
(const mxArray *mxclass, int N, char *fieldname){
    // copy the field to a mxArray
    mxArray *mxField = mxGetProperty(mxclass, 0, fieldname);
    mxArray *mxDouble;
    mexCallMATLAB(1, &mxDouble, 1, &mxField, "int8");
    mxDestroyArray(mxField);
    
    // get the data from the mxArray
    signed char *type = (signed char *)mxGetData(mxDouble);
    NdgCellType *typeVec = (NdgCellType *)calloc(N, sizeof(NdgCellType));
    for (int n=0; n < N; n++)
        typeVec[n] = (NdgCellType)type[n];
    
    mxDestroyArray(mxDouble);
    return typeVec;
}

/** @brief Copy the NdgCell form Matlab variable to C structure. */
NdgCell *mxGetNdgCell(const mxArray *mxCell){
    
    NdgCell *cell = mallocNdgCell();
    // field:: N
    cell->N = (int)mxGetPropertyDouble(mxCell, "N");
    
    // field:: type
    cell->type = (NdgCellType)mxGetPropertyDouble(mxCell, "type");
    // field:: Nv
    cell->Nv = (int)mxGetPropertyDouble(mxCell, "Nv");
    // field:: vr, vs, vt
    cell->vr = mxGetPropertyDoubleVector(mxCell, cell->Nv, "vr");
    cell->vs = mxGetPropertyDoubleVector(mxCell, cell->Nv, "vs");
    cell->vt = mxGetPropertyDoubleVector(mxCell, cell->Nv, "vt");
    // field:: vol
    cell->vol = mxGetPropertyDouble(mxCell, "vol");
    // field:: Nface
    cell->Nface = (int)mxGetPropertyDouble(mxCell, "Nface");
    // field:: Nfv
    cell->Nfv = mxGetPropertyIntVector(mxCell, cell->Nface, "Nfv");
    
    // field:: FToV
    int TFV = 0, maxNFv = 0;
    for (int n = 0; n < cell->Nface; n++)
    {
        TFV += cell->Nfv[n];
        maxNFv = max(maxNFv, cell->Nfv[n]);
    }
    
    cell->FToV = (int **)calloc(cell->Nface, sizeof(int *));
    cell->FToV[0] = (int *)calloc(TFV, sizeof(int));
    for (int n = 0; n < (cell->Nface - 1); n++)
        cell->FToV[n + 1] = cell->FToV[n] + cell->Nfv[n];
    
    int *iTemp = mxGetPropertyIntVector(mxCell, maxNFv * cell->Nface, "FToV");
    for (int f = 0; f < cell->Nface; f++)
        for (int n = 0; n < cell->Nfv[f]; n++)
        {
            cell->FToV[f][n] = iTemp[f * maxNFv + n];
        }
    free(iTemp);
    // field:: faceType
    cell->faceType = mxGetPropertyFaceTypeVector(mxCell, cell->Nface, "faceType");
    
    // field:: Np
    cell->Np = (int)mxGetPropertyDouble(mxCell, "Np");
    
    // field:: r, s, t
    cell->r = mxGetPropertyDoubleVector(mxCell, cell->Np, "r");
    cell->s = mxGetPropertyDoubleVector(mxCell, cell->Np, "s");
    cell->t = mxGetPropertyDoubleVector(mxCell, cell->Np, "t");
    
    // field:: matrix
    cell->M = mxGetPropertyDoubleMatrix(mxCell, cell->Np, cell->Np, "M");
    cell->V = mxGetPropertyDoubleMatrix(mxCell, cell->Np, cell->Np, "V");
    cell->Dr = mxGetPropertyDoubleMatrix(mxCell, cell->Np, cell->Np, "Dr");
    cell->Ds = mxGetPropertyDoubleMatrix(mxCell, cell->Np, cell->Np, "Ds");
    cell->Dt = mxGetPropertyDoubleMatrix(mxCell, cell->Np, cell->Np, "Dt");
    
    cell->Drt = transMatrix(cell->Np, cell->Np, cell->Dr);
    cell->Dst = transMatrix(cell->Np, cell->Np, cell->Ds);
    cell->Dtt = transMatrix(cell->Np, cell->Np, cell->Dt);
    
    // field:: Nq
    cell->Nq = (int)mxGetPropertyDouble(mxCell, "Nq");
    
    // field:: rq, sq, tq
    cell->rq = mxGetPropertyDoubleVector(mxCell, cell->Nq, "rq");
    cell->sq = mxGetPropertyDoubleVector(mxCell, cell->Nq, "sq");
    cell->tq = mxGetPropertyDoubleVector(mxCell, cell->Nq, "tq");
    cell->wq = mxGetPropertyDoubleVector(mxCell, cell->Nq, "tq");
    
    // field:: Vq
    cell->Vq = mxGetPropertyDoubleMatrix(mxCell, cell->Np, cell->Nq, "Vq");
    // field:: Nfp
    cell->Nfp = mxGetPropertyIntVector(mxCell, cell->Nface, "Nfp");
    // field:: TNfp
    cell->TNfp = (int)mxGetPropertyDouble(mxCell, "TNfp");
    
    // field:: Fmask
    int maxNfp = 0;
    for (int n = 0; n < cell->Nface; n++)
        maxNfp = max(maxNfp, cell->Nfp[n]);
    
    cell->Fmask = (int **)calloc(cell->Nface, sizeof(int *));
    cell->Fmask[0] = (int *)calloc(cell->TNfp, sizeof(int));
    for (int n = 0; n < (cell->Nface - 1); n++)
        cell->Fmask[n + 1] = cell->Fmask[n] + cell->Nfp[n];
    
    iTemp = mxGetPropertyIntVector(mxCell, maxNfp * cell->Nface, "Fmask");
    for (int f = 0; f < cell->Nface; f++)
        for (int n = 0; n < cell->Nfp[f]; n++)
            cell->Fmask[f][n] = iTemp[f * maxNfp + n];
    free(iTemp);
    // field:: LIFT
    cell->LIFT = mxGetPropertyDoubleMatrix(mxCell, cell->TNfp, cell->Np, "LIFT");
    cell->LIFTt = transMatrix(cell->TNfp, cell->Np, cell->LIFT);
    return cell;
}

void mxPrintNdgCellInfo(NdgCell *cell, char *message)
{
    printf("%s: %s\n", __FUNCTION__, message);
    mexPrintf("{\n");
    mexPrintf("\t N = %d\n", cell->N);
    mexPrintf("\t type = %d\n", (int)cell->type);
    mexPrintf("\t Nv = %d\n", cell->Nv);
    mexPrintf("\t [ vr, vs, vt ] = [\n");
    for (int n = 0; n < cell->Nv; n++)
        mexPrintf("\t\t %f, %f, %f\n", cell->vr[n], cell->vs[n], cell->vt[n]);
    mexPrintf("\t ]\n");
    
    mexPrintf("\t vol = %f\n", cell->vol);
    mexPrintf("\t Nface = %d\n", cell->Nface);
    mexPrintf("\t Nfv = [");
    for (int n = 0; n < cell->Nface; n++)
        mexPrintf("%d, ", cell->Nfv[n]);
    mexPrintf("]\n");
    
    mexPrintf("\t FToV = [\n");
    for (int f = 0; f < cell->Nface; f++)
    {
        mexPrintf("\t\t ");
        for (int n = 0; n < cell->Nfv[f]; n++)
            mexPrintf("%d, ", cell->FToV[f][n]);
        mexPrintf("\n");
    }
    
    mexPrintf("\t faceType = [");
    for (int n = 0; n < cell->Nface; n++)
        mexPrintf("%d, ", (int)cell->faceType[n]);
    mexPrintf("]\n");
    
    mexPrintf("\t Np = %d\n", cell->Np);
    mexPrintf("\t (r, s, t) = [\n");
    for (int n = 0; n < cell->Np; n++)
        mexPrintf("\t\t %f, %f, %f\n", cell->r[n], cell->s[n], cell->t[n]);
    
    mexPrintDoubleMatrix("\t V = [", cell->V, cell->Np, cell->Np);
    mexPrintf("\t\t]\n");
    
    mexPrintDoubleMatrix("\t M = [", cell->M, cell->Np, cell->Np);
    mexPrintf("\t\t]\n");
    
    mexPrintDoubleMatrix("\t Dr = [", cell->Dr, cell->Np, cell->Np);
    mexPrintf("\t\t]\n");
    
    mexPrintDoubleMatrix("\t Ds = [", cell->Ds, cell->Np, cell->Np);
    mexPrintf("\t\t]\n");
    
    mexPrintDoubleMatrix("\t Drt = [", cell->Drt, cell->Np, cell->Np);
    mexPrintf("\t\t]\n");
    
    mexPrintDoubleMatrix("\t Dst = [", cell->Dst, cell->Np, cell->Np);
    mexPrintf("\t\t]\n");
    
    mexPrintf("\t Nq = %d\n", cell->Nq);
    mexPrintf("\t [ rq, sq, tq, wq ] = [\n");
    for (int n = 0; n < cell->Nq; n++)
        mexPrintf("\t\t %f, %f, %f, %f\n",
                  cell->rq[n], cell->sq[n], cell->tq[n], cell->wq[n]);
    mexPrintf("\t ]\n");
    
    mexPrintDoubleMatrix("\t Vq = [", cell->Vq, cell->Np, cell->Nq);
    mexPrintf("\t\t]\n");
    
    mexPrintf("\t Nfp = [");
    for (int n = 0; n < cell->Nface; n++)
        mexPrintf("%d, ", cell->Nfp[n]);
    mexPrintf("]\n");
    mexPrintf("\t TNfp = %d\n", cell->TNfp);
    
    mexPrintf("\t Fmask = [\n");
    for (int f = 0; f < cell->Nface; f++)
    {
        mexPrintf("\t\t ");
        for (int n = 0; n < cell->Nfp[f]; n++)
            mexPrintf("%d, ", cell->Fmask[f][n]);
        mexPrintf("\n");
    }
    
    mexPrintDoubleMatrix("\t LIFT = [", cell->LIFT, cell->TNfp, cell->Np);
    mexPrintf("\t\t]\n");
    mexPrintDoubleMatrix("\t LIFTt = [", cell->LIFTt, cell->Np, cell->TNfp);
    mexPrintf("\t\t]\n");
    
    mexPrintf("}\n");
    return;
}