//
//  AdvSolver2d.c
//  NDG-FEM
//
//  Created by li12242 on 17/10/17.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#include "NdgSolver.h"
#include <time.h>

typedef struct AdvectionSolver2d{
    NdgSolver *solver;
    double **u;
    double **v;
} AdvectionSolver2d;

AdvectionSolver2d adv;

void evaluateAdvectionNumericalFlux(int ind,
                                    double *fM, double *fP,
                                    double nx, double ny, double nz,
                                    double *flux){
    const double unM = adv.u[0][ind] * nx + adv.v[0][ind] * ny;
    if (unM > 0){
        *flux = fM[0] * unM;
    }else{
        *flux = fP[0] * unM;
    }
    return;
}

void evaluateAdvectionNodalFlux(int ind,
                                double *f_Q,
                                double *E, double *G, double *H){
    E[0] = f_Q[0] * adv.u[0][ind];
    G[0] = f_Q[0] * adv.v[0][ind];
    return;
}

void evaluateAdvectionNeighbourValue(double *fM, double *fP, double *fext,
                                     double nx, double ny, double nz,
                                     NdgEdgeType edgeType,
                                     double *f_extP){
    switch (edgeType)
    {
    case NdgEdgeInner:
        f_extP[0] = fP[0];
        break;
    default:
        f_extP[0] = 0.;
        break;
    }
    return;
}

void evaluateAdvectionSourceTerm(NdgPhysUnion *phys, double ***var, double ***rhs){
    return;
}

void evaluatePostFunc(NdgPhysUnion *phys, double ***rhs){
    return;
}

void outputPhysResult(NdgOutputFileInfo *outputFileInfo, NdgPhysUnion *phys, double time, int tstep){
    NdgNcFile *ncfile = outputFileInfo->ncfile;
    size_t start_t = tstep;
    size_t count_t = 1;
    size_t start_f[4] = {tstep, 0, 0, 0};
    size_t count_f[4] = {1, phys->Nfield, phys->mesh->K, phys->mesh->cell->Np};
    // output time
    nc_put_vara_double(ncfile->ncid, ncfile->var[0].id, &start_t, &count_t, &time);
    // output physical field
    nc_put_vara_double(ncfile->ncid, ncfile->var[1].id, start_f, count_f, phys->fvar[0][0]);
    return;
}

static void mxCopyResult(mxArray *mxfield){
    NdgPhysUnion *phys = adv.solver->phys;
    const int Nfield = phys->Nfield;
    const int K = phys->mesh->K;
    const int Np = phys->mesh->cell->Np;

    double *field_ptr = mxGetPr(mxfield);
    int sk = 0;
    for (int fld = 0; fld < Nfield; fld++){
        for (int k = 0; k < K; k++){
            for (int n = 0; n < Np; n++){
                field_ptr[sk++] = phys->fvar[fld][k][n];
            }
        }
    }
    return;
}

static void freeSolver(){
    freeNdgSolver(adv.solver);
    freeMatrix(adv.u);
    freeMatrix(adv.v);
    return;
}

void mxGetSolver(const mxArray *mxSolver){
    mxArray *mxPhys = mxGetProperty(mxSolver, 0, "phys");
    NdgPhysUnion *phys = mxGetNdgPhysUnion(mxPhys);
    mxDestroyArray(mxPhys);
    mxCopyInitialPhysUnionField(mxSolver, "fvar", phys);
    const int K = phys->mesh->K;
    const int Np = phys->mesh->cell->Np;
    adv.u = mxGetPropertyDoubleMatrix(mxSolver, K, Np, "u");
    adv.v = mxGetPropertyDoubleMatrix(mxSolver, K, Np, "v");
    
    adv.solver = initNdgSolver
    (phys,
     NULL,
     NULL,
     outputPhysResult,
     NdgRHS2d,
     evaluateAdvectionNodalFlux,
     evaluateAdvectionNumericalFlux,
     evaluateAdvectionNeighbourValue,
     evaluateAdvectionSourceTerm,
     evaluatePostFunc);
    
    return;
}

void mexFunction(int nlhs,
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[]){

    mxGetSolver(prhs[0]);
    
#ifdef PROFILE
    clock_t tstart = clock()/CLOCKS_PER_SEC;
#endif
    
    adv.solver->evaluateTemporalDiscrete(adv.solver);

#ifdef PROFILE
    clock_t tend = clock()/CLOCKS_PER_SEC;
    
    extern double elapsedTimeSurfaceTerm;
    extern double elapsedTimeVolumeTerm;
    
    mexPrintf("Total elapsed time: %f\n", (double)(tend-tstart) );
    mexPrintf("Surface term elapsed time: %f\n", (double)elapsedTimeSurfaceTerm );
    mexPrintf("Volume term elapsed time: %f\n", (double)elapsedTimeVolumeTerm );
#endif
    
    plhs[0] = mxGetProperty(prhs[0], 0, "fvar");
    mxCopyResult(plhs[0]);
    freeSolver();
    return;
}