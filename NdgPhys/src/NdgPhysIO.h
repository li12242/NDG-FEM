//
//  NdgPhysIO.h
//  NDG-FEM
//
//  Created by li12242 on 17/10/19.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#ifndef __NDG_FEM__NdgPhysIO__
#define __NDG_FEM__NdgPhysIO__

#include "NdgPhys.h"

typedef struct NdgOutputFileInfo{
    NdgNcFile *ncfile;
    int ioStep;
    double previousOutputTime;
    int previousOutputStep;
    void (*outputPhysResult)(struct NdgOutputFileInfo *outputFileInfo, NdgPhysUnion *phys, double time, int tstep);
}NdgOutputFileInfo;

/** */
void updateExternalInfoFromFile(NdgPhysUnion *phys, double time, int tstep);

/** */
void updateOutputResultToFile();

/** */
void updateOutputResultByDeltaTime
(NdgPhysUnion *phys, NdgOutputFileInfo *outputFileInfo, double time, int tstep);

/** */
void updateOutputResultByDeltaStep
(NdgPhysUnion *phys, NdgOutputFileInfo *outputFileInfo, double time, int tstep);

/** */
NdgOutputFileInfo* initOutputFileInfo
(NdgPhysUnion *phys,
 void (*updatePhysResult)(NdgOutputFileInfo* outputFileInfo, NdgPhysUnion *phys, double time, int tstep) );

/** */
void freeOutputFileInfo(NdgOutputFileInfo* outputFileInfo);

#endif /* defined(__NDG_FEM__NdgPhysIO__) */
