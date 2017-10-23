//
//  UMeshIO.c
//  NDGOM_Mex
//
//  Created by li12242 on 17/10/10.
//  Copyright (c) 2017年 li12242. All rights reserved.
//
#include "NdgMeshIO.h"

/** Allocate the memory for NdgMeshUnion structure. */
static NdgUMeshUnion *mallocNdgUMeshUnion(){
    NdgUMeshUnion *mesh = (NdgUMeshUnion *)calloc(1, sizeof(NdgUMeshUnion));
    return mesh;
}

/** Get vector of NdgRegionType from Matlab class variable. */
static NdgRegionType *mxGetPropertyRegionTypeVector(const mxArray *mxclass, int N, char *fieldname)
{
    mxArray *mxField = mxGetProperty(mxclass, 0, fieldname);
    mxArray *mxDouble;
    mexCallMATLAB(1, &mxDouble, 1, &mxField, "int8");
    mxDestroyArray(mxField);
    signed char *type;
    type = (signed char *)mxGetData(mxDouble);
    NdgRegionType *typeVec = (NdgRegionType *)calloc(N, sizeof(NdgRegionType));
    for (int n = 0; n < N; n++)
        typeVec[n] = (NdgRegionType)type[n];
    
    mxDestroyArray(mxDouble);
    return typeVec;
}

static
NdgEdgeType**
mxGetPropertyEdgeTypeMatrix(const mxArray *mxclass,
                            int Nrow,
                            int Ncol,
                            char *fieldname)
{
    mxArray *mxField = mxGetProperty(mxclass, 0, fieldname);
    mxArray *mxDouble;
    mexCallMATLAB(1, &mxDouble, 1, &mxField, "int8");
    mxDestroyArray(mxField);
    signed char *type = (signed char *)mxGetData(mxDouble);
    NdgEdgeType **typeVec = (NdgEdgeType **)calloc( Nrow, sizeof(NdgEdgeType*) );
    typeVec[0] = (NdgEdgeType *)calloc( Nrow*Ncol, sizeof(NdgEdgeType) );
    for (int n = 1; n < Nrow; ++n)
        typeVec[n] = typeVec[n - 1] + Ncol;
    
    for (int n=0; n<Nrow*Ncol; n++)
        typeVec[0][n] = (NdgEdgeType)type[n];
    
    mxDestroyArray(mxDouble);
    return typeVec;
}

/** Copy the NdgMeshUnion structure from Matlab variable. */
NdgUMeshUnion *
mxGetNdgUMeshUnion(const mxArray *mxMesh){
    
    NdgUMeshUnion *mesh = mallocNdgUMeshUnion();
    mxArray *mxField = mxGetProperty(mxMesh, 0, "cell");
    mesh->cell = mxGetNdgCell(mxField);
    mxDestroyArray(mxField);
    
    mxField = mxGetProperty(mxMesh, 0, "uedge");
    mxDestroyArray(mxField);
    
    mesh->K = (int)mxGetPropertyDouble(mxMesh, "K");
    mesh->Nv = (int)mxGetPropertyDouble(mxMesh, "Nv");
    mesh->EToV = mxGetPropertyIntMatrix(mxMesh, mesh->K, mesh->cell->Nv, "EToV");
    mesh->EToE = mxGetPropertyIntMatrix(mxMesh, mesh->K, mesh->cell->Nface, "EToE");
    mesh->EToF = mxGetPropertyIntMatrix(mxMesh, mesh->K, mesh->cell->Nface, "EToF");
    mesh->EToM = mxGetPropertyIntMatrix(mxMesh, mesh->K, mesh->cell->Nface, "EToM");
    mesh->EToR = mxGetPropertyRegionTypeVector(mxMesh, mesh->K, "EToR");
    mesh->EToB = mxGetPropertyEdgeTypeMatrix(mxMesh, mesh->K, mesh->cell->Nface, "EToB");
    mesh->vx = mxGetPropertyDoubleVector(mxMesh, mesh->Nv, "vx");
    mesh->vy = mxGetPropertyDoubleVector(mxMesh, mesh->Nv, "vy");
    mesh->vz = mxGetPropertyDoubleVector(mxMesh, mesh->Nv, "vz");
    
    mesh->x = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "x");
    mesh->y = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "y");
    mesh->z = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "z");
    mesh->J = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "J");
    
    mesh->rx = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "rx");
    mesh->ry = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "ry");
    mesh->rz = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "rz");
    mesh->sx = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "sx");
    mesh->sy = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "sy");
    mesh->sz = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "sz");
    mesh->tx = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "tx");
    mesh->ty = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "ty");
    mesh->tz = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->Np, "tz");
    
    mesh->eidM = mxGetPropertyIntMatrix(mxMesh, mesh->K, mesh->cell->TNfp, "eidM");
    mesh->eidP = mxGetPropertyIntMatrix(mxMesh, mesh->K, mesh->cell->TNfp, "eidP");
    mesh->eidType = mxGetPropertyEdgeTypeMatrix(mxMesh, mesh->K, mesh->cell->Nface, "eidtype");
    mesh->Js = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->TNfp, "Js");
    mesh->nx = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->TNfp, "nx");
    mesh->ny = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->TNfp, "ny");
    mesh->nz = mxGetPropertyDoubleMatrix(mxMesh, mesh->K, mesh->cell->TNfp, "nz");
    return mesh;
}

/**  */
void
mxPrintUMeshUnion(NdgUMeshUnion *mesh, char *message){
    
    printf("%s: %s\n", __FUNCTION__, message);
    mexPrintf("{\n");
    mexPrintf("\t K = %d\n", mesh->K);
    mexPrintf("\t Nv = %d\n", mesh->Nv);
    mexPrintIntMatrix("\t EToV = [", mesh->EToV, mesh->K, mesh->cell->Nv);
    mexPrintf("\t ]");
    mexPrintIntMatrix("\t EToE = [", mesh->EToE, mesh->K, mesh->cell->Nv);
    mexPrintf("\t ]");
    mexPrintIntMatrix("\t EToF = [", mesh->EToF, mesh->K, mesh->cell->Nv);
    mexPrintf("\t ]");
    
    mexPrintf("\t EToR = [\n\t\t");
    for (int n = 0; n < mesh->K; n++)
        mexPrintf("%d, ", (int)mesh->EToR[n]);
    mexPrintf("\t ]\n");
    
    mexPrintf("\t [vx, vy, vz] = [\t");
    for (int n = 0; n < mesh->Nv; n++)
        mexPrintf("\t %f, %f, %f\n", mesh->vx[n], mesh->vy[n], mesh->vz[n]);
    mexPrintf("\t ]\n");
    
    mexPrintDoubleMatrix("\t x = [", mesh->x, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t y = [", mesh->y, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t z = [", mesh->z, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t J = [", mesh->J, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t rx = [", mesh->rx, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t ry = [", mesh->ry, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t rz = [", mesh->rz, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t sx = [", mesh->sx, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t sy = [", mesh->sy, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t sz = [", mesh->sz, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t tx = [", mesh->tx, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t ty = [", mesh->ty, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t tz = [", mesh->tz, mesh->K, mesh->cell->Np);
    mexPrintf("\t ]");
    
    mexPrintf("\t EToB = [\n\t\t");
    for (int n = 0; n < mesh->K; n++)
        for (int m = 0; m < mesh->cell->Nface; m++) {
            mexPrintf("%d, ", (int)mesh->EToB[n][m]);
        }
    
    mexPrintf("\t eidtype = [\n\t\t");
    for (int n = 0; n < mesh->K; n++)
        for (int m = 0; m < mesh->cell->TNfp; m++) {
            mexPrintf("%d, ", (int)mesh->eidType[n][m]);
        }
    mexPrintf("\t ]\n");
    mexPrintIntMatrix("\t eidM = [", mesh->eidM, mesh->K, mesh->cell->TNfp);
    mexPrintf("\t ]");
    mexPrintIntMatrix("\t eidP = [", mesh->eidP, mesh->K, mesh->cell->TNfp);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t nx = [", mesh->nx, mesh->K, mesh->cell->TNfp);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t ny = [", mesh->ny, mesh->K, mesh->cell->TNfp);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t nz = [", mesh->nz, mesh->K, mesh->cell->TNfp);
    mexPrintf("\t ]");
    mexPrintDoubleMatrix("\t Js = [", mesh->Js, mesh->K, mesh->cell->TNfp);
    mexPrintf("\t ]");
    mexPrintf("}\n");
}

void freeUMeshUnionOutpufFile(NdgNcFile *ncfile){
    free(ncfile->dim);
    free(ncfile->var);
    free(ncfile);
}

void freeBasicUMeshUnionOutputFile(NdgNcFile *ncfile){
    free(ncfile->dim);
    free(ncfile->var);
    closeNetcdfFile(ncfile[0]);
    free(ncfile);
}

NdgNcFile* setupBasicUMeshUnionOutputFile(NdgUMeshUnion mesh, char *casename, int Nfield){
    
    NdgCell *cell = mesh.cell;
    const int N = cell->N;
    const int Np = cell->Np;
    const int Nface = cell->Nface;
    const int TNfp = cell->TNfp;
    const int K = mesh.K;
    /* define dimensions */
    const int NDim = 11;
    NdgNcDim *ncdim = (NdgNcDim*) calloc(NDim, sizeof(NdgNcDim));
    
    setNcDim(ncdim    , "N" , N);
    setNcDim(ncdim + 1, "Np", Np);
    setNcDim(ncdim + 2, "CellType", cell->type);
    setNcDim(ncdim + 3, "Nq", cell->Nq);
    setNcDim(ncdim + 4, "Nv", cell->Nv);
    setNcDim(ncdim + 5, "Nface", cell->Nface);
    setNcDim(ncdim + 6, "TNfp", TNfp);
    setNcDim(ncdim + 7, "K", K);
    setNcDim(ncdim + 8, "Nvert", mesh.Nv);
    setNcDim(ncdim + 9, "Nfield", Nfield);
    setNcDim(ncdim + 10, "Time", 0);

    /* define variables */
    const int NVar = 8;
    NdgNcVar *ncvar = (NdgNcVar*) calloc(NVar, sizeof(NdgNcVar));
    
    int fieldDimId[4] = {10, 9, 7, 1};
    setNcVar(ncvar    , "time", 1, fieldDimId, NC_DOUBLE);
    setNcVar(ncvar + 1, "fvar", 4, fieldDimId, NC_DOUBLE);
    int meshDimId[2] = {7, 4};
    setNcVar(ncvar + 2, "EToV", 2, meshDimId, NC_INT);
    setNcVar(ncvar + 3, "EToR", 1, meshDimId, NC_INT);
    int bcDimId[2] = {7, 5};
    setNcVar(ncvar + 4, "EToB", 2, bcDimId, NC_INT);
    int vertDimId[1] = {8};
    setNcVar(ncvar + 5, "vx", 1, vertDimId, NC_DOUBLE);
    setNcVar(ncvar + 6, "vy", 1, vertDimId, NC_DOUBLE);
    setNcVar(ncvar + 7, "vz", 1, vertDimId, NC_DOUBLE);


    NdgNcFile *ncfile = (NdgNcFile *)calloc(1, sizeof(NdgNcFile));
    setNcFile(ncfile, casename, NDim, ncdim, NVar, ncvar);
    /* generate output NetCDF file */
    makeNetcdfFile(ncfile);
    
    /* put mesh variables */
    int varId = 2;
    nc_put_var_int(ncfile->ncid, ncvar[varId++].id, mesh.EToV[0]);
    
    int EToR[K];
    for (int k = 0; k<K; k++ ) EToR[k] = (int)mesh.EToR[k];
    nc_put_var_int(ncfile->ncid, ncvar[varId++].id, EToR);
    
    int EToB[K*Nface], sk=0;
    for (int k=0; k<K; k++ )
        for (int f=0; f<Nface; f++) {
            EToB[sk++] = (int)mesh.EToB[k][f];
        }
    nc_put_var_int(ncfile->ncid, ncvar[varId++].id, EToB);
    nc_put_var_double(ncfile->ncid, ncvar[varId++].id, mesh.vx);
    nc_put_var_double(ncfile->ncid, ncvar[varId++].id, mesh.vy);
    nc_put_var_double(ncfile->ncid, ncvar[varId++].id, mesh.vz);
    
    return ncfile;
}

void outputBasicUMeshUnionResult(NdgNcFile ncfile, int tstep, double time, double ***fvar){
    const int K = ncfile.dim[7].length;
    const int Np = ncfile.dim[1].length;
    const int Nfield = ncfile.dim[9].length;
    size_t start_v[4] = {tstep, 0, 0, 0};
    size_t count_v[4] = {1, Nfield, K, Np};
    size_t start_t = tstep;
    size_t count_t = 1;
    
    nc_put_vara_double(ncfile.ncid, ncfile.var[0].id, &start_t, &count_t, &time);
    nc_put_vara_double(ncfile.ncid, ncfile.var[1].id, start_v, count_v, fvar[0][0]);
    
    return;
}