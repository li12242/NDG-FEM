//
//  NdgPhysIO.c
//  NDG-FEM
//
//  Created by li12242 on 17/10/19.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#include "NdgPhysIO.h"

void updateOutputResultByDeltaTime
(NdgPhysUnion *phys, NdgOutputFileInfo *outputFileInfo, double time, int tstep)
{
    if (( time -  outputFileInfo->previousOutputTime) > phys->outputTimeInterval){
        outputFileInfo->outputPhysResult(outputFileInfo, phys, time, outputFileInfo->ioStep);
        outputFileInfo->ioStep++;
        outputFileInfo->previousOutputTime = time;
    }
    return;
}

void updateOutputResultByDeltaStep
(NdgPhysUnion *phys, NdgOutputFileInfo *outputFileInfo, double time, int tstep)
{
    if (( tstep -  outputFileInfo->previousOutputStep) > phys->outputStepInterval){
        outputFileInfo->outputPhysResult(outputFileInfo, phys, time, outputFileInfo->ioStep);
        outputFileInfo->ioStep++;
        outputFileInfo->previousOutputStep = tstep;
    }
    return;
}

void freeOutputFileInfo(NdgOutputFileInfo* outputFileInfo){
    freeBasicUMeshUnionOutputFile(outputFileInfo->ncfile);
    free(outputFileInfo);
}

NdgOutputFileInfo* initOutputFileInfo
(NdgPhysUnion *phys,
 void (*outputPhysResult)(NdgOutputFileInfo* outputFileInfo, NdgPhysUnion *phys, double time, int tstep) )
{
    NdgOutputFileInfo *outputFileInfo = (NdgOutputFileInfo *)calloc(1, sizeof(NdgOutputFileInfo));
    // read the output
    outputFileInfo->ncfile = setupBasicUMeshUnionOutputFile(phys->mesh[0], phys->outputNetcdfCaseName, phys->Nfield);
    outputFileInfo->ioStep = 0;
    outputFileInfo->previousOutputStep = 0;
    outputFileInfo->previousOutputTime = 0.0;
    outputFileInfo->outputPhysResult = outputPhysResult;
    return outputFileInfo;
}

/** update external value from NetCDF file */
void updateExternalInfoFromFile (NdgPhysUnion *phys, double time, int tstep){
    return;
}