//
//  mxNdgPhys.c
//  NDG-FEM
//
//  Created by li12242 on 17/10/21.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#include "mxNdgPhys.h"

static NdgPhysUnion *mallocNdgPhysUnion(){
    NdgPhysUnion *phys = (NdgPhysUnion *)calloc(1, sizeof(NdgPhysUnion));
    return phys;
}

void
mxCopyInitialPhysUnionField
(const mxArray *mxPhys,
 const char *fieldName,
 NdgPhysUnion *phys)
{
    mxArray *mxField = mxGetProperty(mxPhys, 0, fieldName);
    double *field_ptr = mxGetPr(mxField);
    
    double ***fvar = phys->fvar;
    const int Nfield = phys->Nfield;
    const int K = phys->mesh->K;
    const int Np = phys->mesh->cell->Np;

    int ind = 0;
    for (int m = 0; m < Nfield; m++){
        for (int k = 0; k<K; k++) {
            for (int n = 0; n<Np; n++) {
                fvar[m][k][n] = field_ptr[ind++];
            }
        }
    }
    mxDestroyArray(mxField);
}

static void mxGetlimiterInfo(const mxArray *mxPhys, NdgPhysUnion *phys){
    double type = mxGetPropertyDouble(mxPhys, "limiterType");
    phys->limiterType = (NdgLimiterType)type;
    return;
}

static void mxGetTemporalDiscreteInfo(const mxArray *mxPhys, NdgPhysUnion *phys){
    double type = mxGetPropertyDouble(mxPhys, "temporalDiscreteType");
    phys->temporalDiscreteType = (NdgTemporalDiscreteType) type;
    
    type = mxGetPropertyDouble(mxPhys, "timeIntervalType");
    phys->timeIntervalType = (NdgIntervalType) type;
    if ( phys->timeIntervalType == NdgIntervalConst ){
        phys->timeInterval = mxGetPropertyDouble(mxPhys, "timeInterval");
    }
    return;
}

static void mxGetOBCInfo(const mxArray *mxPhys, NdgPhysUnion *phys){
    
    double type = mxGetPropertyDouble(mxPhys, "obcType");
    phys->obcType = (NdgBCType) type;
    if (phys->obcType != NdgBCNone) {
        
        type = mxGetPropertyDouble(mxPhys, "obcIntervalType");
        phys->obcIntervalType = (NdgIntervalType) type;
        if (phys->obcIntervalType == NdgIntervalDeltaTime) {
            phys->obcTimeInterval = mxGetPropertyDouble(mxPhys, "obcTimeInterval");
        }else if(phys->obcIntervalType == NdgIntervalDeltaStep) {
            phys->obcStepInterval = (int) mxGetPropertyDouble(mxPhys, "obcStepInterval");
        }
    }

    mxArray *mxfield = mxGetProperty(mxPhys, 0, "obcFileName");
    mxGetString(mxfield, phys->obcFileName, NC_MAX_NAME_LEN);
    mxDestroyArray(mxfield);
    return;
}

static void mxGetOutputInfo(const mxArray *mxPhys, NdgPhysUnion *phys){
    mxArray *mxfield = mxGetProperty(mxPhys, 0, "outputNetcdfFileName");
    mxGetString(mxfield, phys->outputNetcdfCaseName, NC_MAX_NAME_LEN);
    mxDestroyArray(mxfield);
    
    mxfield = mxGetProperty(mxPhys, 0, "outputNetcdfFileName");
    mxGetString(mxfield, phys->outputNetcdfCaseName, NC_MAX_NAME_LEN);
    mxDestroyArray(mxfield);
    
    double type = mxGetPropertyDouble(mxPhys, "outputIntervalType");
    phys->outputIntervalType = (NdgIntervalType)type;
    if (phys->outputIntervalType == NdgIntervalDeltaTime) {
        phys->outputTimeInterval = mxGetPropertyDouble(mxPhys, "outputInterval");
    }else if (phys->outputIntervalType == NdgIntervalDeltaStep) {
        phys->outputStepInterval = (int) mxGetPropertyDouble(mxPhys, "outputInterval");
    }
    return;
}

static void mxGetMesh(const mxArray *mxPhys, NdgPhysUnion *phys){
    mxArray *mxfield = mxGetProperty(mxPhys, 0, "mesh");
    phys->mesh = mxGetNdgUMeshUnion(mxfield);
    mxDestroyArray(mxfield);
    return;
}

NdgPhysUnion *mxGetNdgPhysUnion(const mxArray *mxPhys){
    NdgPhysUnion *phys = mallocNdgPhysUnion();
    
    // field: limiterType
    mxGetlimiterInfo(mxPhys, phys);
    // field: temporal discrete
    mxGetTemporalDiscreteInfo(mxPhys, phys);
    
    // field: obc type
    mxGetOBCInfo(mxPhys, phys);
    
    // field:: output
    mxGetOutputInfo(mxPhys, phys);
    
    // field:: Nfield & finalTime
    phys->Nfield = (int) mxGetPropertyDouble(mxPhys, "Nfield");
    phys->finalTime = mxGetPropertyDouble(mxPhys, "finalTime");
    // filed:: mesh
    mxGetMesh(mxPhys, phys);
    // field:: variables
    NdgUMeshUnion *mesh = phys->mesh;
    phys->fvar = mallocMatrixVector(phys->Nfield, mesh->K, mesh->cell->Np);
    phys->fext = mallocMatrixVector(phys->Nfield, mesh->K, mesh->cell->Np);
    phys->frhs = mallocMatrixVector(phys->Nfield, mesh->K, mesh->cell->Np);
    return phys;
}

void mxPrintNdgPhys(NdgPhysUnion *phys, char *message){
    mexPrintf("%s: %s\n", __FUNCTION__, message);
    mexPrintf("{\n");
    mexPrintf("\t Nfield = %d\n", phys->Nfield);
    mexPrintf("\t finalTime = %f\n", phys->finalTime);
    mexPrintf("\t limiterType = %d\n", (int) phys->limiterType);
    
    mexPrintf("\t temporalDiscreteType = %d\n", (int)phys->temporalDiscreteType);
    mexPrintf("\t timeIntervalType = %d\n", (int)phys->timeIntervalType);
    mexPrintf("\t timeInterval = %f\n", phys->timeInterval);
    
    mexPrintf("\t obcType = %d\n", (int)phys->obcType);
    mexPrintf("\t obcIntervalType = %d\n", (int)phys->obcIntervalType);
    mexPrintf("\t obcTimeInterval = %f\n", phys->obcTimeInterval);
    mexPrintf("\t obcStepInterval = %d\n", phys->obcStepInterval);
    mexPrintf("\t obcFileName = %s\n", phys->obcFileName);
    
    mexPrintf("\t outputNetcdfFileName = %s\n", phys->outputNetcdfCaseName);
    mexPrintf("\t outputIntervalType = %d\n", (int)phys->outputIntervalType);
    mexPrintf("\t outputStepInterval = %d\n", phys->outputStepInterval);
    mexPrintf("\t outputTimeInterval = %f\n", phys->outputTimeInterval);
    
    for (int m = 0; m<phys->Nfield; m++) {
        mexPrintDoubleMatrix("\t fld = [", phys->fvar[m], phys->mesh->K, phys->mesh->cell->Np);
        mexPrintf("\t ]");
    }
    mexPrintf("}\n");
}