//
//  NdgPhysTest.c
//  NDG-FEM
//
//  Created by li12242 on 17/10/17.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#include "NdgPhys.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
    NdgPhysUnion *phys = mxGetNdgPhysUnion(prhs[0]);
    mxPrintNdgPhys(phys, "Phys");
    freeNdgPhys(phys);
    return;
}