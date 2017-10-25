//
//  NdgPhysTest.c
//  NDG-FEM
//
//  Created by li12242 on 17/10/17.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#include "NdgSolver.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
    NdgSolver solver;
    mxGetNdgSolver(prhs[0], &solver);
    mxPrintNdgSolver(solver, "Solver");
    return;
}