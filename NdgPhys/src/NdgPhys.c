//
//  NdgPhys.c
//  NDG-FEM
//
//  Created by li12242 on 17/10/16.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#include "NdgPhys.h"

void freeNdgPhys(NdgPhysUnion *phys){
    freeNdgUMeshUnion(phys->mesh);
    freeMatrixVector(phys->fvar);
    freeMatrixVector(phys->frhs);
    freeMatrixVector(phys->fext);
    free(phys);
}