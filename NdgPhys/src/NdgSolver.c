//
//  NdgPhysSolver.c
//  NDG-FEM
//
//  Created by li12242 on 17/10/17.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#include "NdgSolver.h"
#include "NdgLimiter.h"
#include "NdgTerporalDiscrete.h"
#include "NdgPhysIO.h"

#ifdef _OPENMP
#include "libiomp/omp.h"
#endif

#ifdef PROFILE
#include <time.h>
double elapsedTimeVolumeTerm = 0.0;
double elapsedTimeSurfaceTerm = 0.0;
#endif

/** return const time interval - dt */
void updatePhysTimeIntervalConstant (NdgPhysUnion *phys, double *dt){
    *dt = phys->timeInterval;
    return;
}

/** */
void updateNoneExternalInfo (NdgPhysUnion *phys, double time, int tstep){
    return;
}

/**  */
void evaluateLimiterNone (NdgPhysUnion *phys, double ***val){
    return;
}

static void initLimiterFunc(NdgSolver *solver, NdgPhysUnion *phys){
    if (phys->limiterType == NdgLimiterNone){
        solver->evaluateLimiter = evaluateLimiterNone;
    }else if (phys->limiterType == NdgLimiterVert){
        solver->evaluateLimiter = vert_limiter;
    }
}

static void initTimeIntervalFunc
(NdgSolver *solver,
 NdgPhysUnion *phys,
 void (*updatePhysTimeInterval)(NdgPhysUnion *phys, double *dt) )
{
    if (phys->timeIntervalType == NdgIntervalConst){
        solver->updatePhysTimeInterval = updatePhysTimeIntervalConstant;
    }else{
        solver->updatePhysTimeInterval = updatePhysTimeInterval;
    }
}

static void initUpdateExternalFileFunc
(NdgSolver *solver,
 NdgPhysUnion *phys,
 void (*updatePhysExternal)(NdgPhysUnion *phys, double time, int tstep))
{
    if (phys->obcType == NdgBCNone) {
        solver->updatePhysExternalInfo = updateNoneExternalInfo;
    }else if(phys->obcType == NdgBCFile){
        solver->updatePhysExternalInfo = updateExternalInfoFromFile;
    }else if(phys->obcType == NdgBCFunc){
        solver->updatePhysExternalInfo = updatePhysExternal;
    }
}

static void initUpdateOutputResultFunc
(NdgSolver *solver,
 NdgPhysUnion *phys,
 void (*outputPhysResult)(NdgOutputFileInfo *outputFileInfo, NdgPhysUnion *phys, double time, int tstep) )
{
    solver->outFileInfo = initOutputFileInfo(phys, outputPhysResult);
    
    if (phys->outputIntervalType == NdgIntervalDeltaTime) {
        solver->updateOutputResult = updateOutputResultByDeltaTime;
    }else if ( phys->outputIntervalType == NdgIntervalDeltaStep ){
        solver->updateOutputResult = updateOutputResultByDeltaStep;
    }
}

NdgSolver* initNdgSolver
(NdgPhysUnion *phys,
 void (*updatePhysTimeInterval)(NdgPhysUnion *phys, double *dt),
 void (*updatePhysExternal)(NdgPhysUnion *phys, double time, int tstep),
 void (*outputPhysResult)(NdgOutputFileInfo *outputFileInfo, NdgPhysUnion *phys, double time, int tstep),
 void (*evaluateRHSTerm)(NdgSolver *solver, NdgPhysUnion *phys, double ***rhs),
 void (*evaluateRHSNodalFlux)(int ind, double *val, double *E, double *G, double *H),
 void (*evaluateRHSNumericalFlux)(int ind, double *fM, double *fP, double nx, double ny, double nz, double *flux),
 void (*evaluateRHSNeighbourValue)(double *fM, double *fP, double *fext, double nx, double ny, double nz, NdgEdgeType edgeNodeType, double *fextP),
 void (*evaluateRHSSourceTerm)(NdgPhysUnion *phys, double ***var, double ***rhs),
 void (*evaluatePostFunc)(NdgPhysUnion *phys, double ***val)
 ){
    
    NdgSolver *solver = (NdgSolver *)calloc(1, sizeof(NdgSolver));
    
    solver->phys = phys;
    
    // set the way to update the computational time interval - dt
    initTimeIntervalFunc(solver, phys, updatePhysTimeInterval);
    
    // set the update_obc_ext function
    initUpdateExternalFileFunc(solver, phys, updatePhysExternal);
    
    // set the limiter
    initLimiterFunc(solver, phys);
    
    // temporal discrete scheme
    if (phys->temporalDiscreteType == NdgTemporalDiscreteRK45) {
        solver->evaluateTemporalDiscrete = SSP_RK45;
    }else if(phys->temporalDiscreteType == NdgTemporalDiscreteEuler){
        solver->evaluateTemporalDiscrete = Euler;
    }
    
    // set the output function
    initUpdateOutputResultFunc(solver, phys, outputPhysResult);
    
    // functions to calculate the rhs term
    solver->evaluateRHSNodalFlux = evaluateRHSNodalFlux;
    solver->evaluateRHSNumericalFlux = evaluateRHSNumericalFlux;
    solver->evaluateRHSNeighbourValue = evaluateRHSNeighbourValue;
    
    solver->evaluateRHSTerm = evaluateRHSTerm;
    solver->evaluateRHSSourceTerm = evaluateRHSSourceTerm;
    solver->evaluatePostFunc = evaluatePostFunc;
    
    return solver;
}

void freeNdgSolver(NdgSolver *solver){
    freeOutputFileInfo(solver->outFileInfo);
    free(solver);
};

static void NdgVolumeTerm2d
(NdgUMeshUnion *mesh,
 void (*evaluateRHSNodalFlux)(int ind, double *val, double *E, double *G, double *H),
 int Nfield,
 double ***f_valQ,
 double ***f_rhsQ){
    
    NdgCell *cell = mesh->cell;
    const int K = mesh->K;
    const int Np = cell->Np;
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k<K; k++) {
        double f_Q[Nfield*Np], E[Nfield*Np], G[Nfield*Np];
        
        for (int n=0; n<Np; n++) {
            int sk = k*Np + n;
            for (int fld=0; fld<Nfield; fld++)
                f_Q[n*Nfield + fld] = f_valQ[fld][k][n];
            
            evaluateRHSNodalFlux(sk, f_Q+n*Nfield, E+n*Nfield, G+n*Nfield, NULL);
        }
        
        for (int n=0; n<Np; n++) {
            double *Dr = cell->Drt[n];
            double *Ds = cell->Dst[n];
            const double drdx = mesh->rx[k][n];
            const double drdy = mesh->ry[k][n];
            const double dsdx = mesh->sx[k][n];
            const double dsdy = mesh->sy[k][n];
            
            double rhs[Nfield];
            for (int fld=0; fld<Nfield; fld++)
                rhs[fld] = 0.0;
            
            int sk = 0;
            for (int m=0; m<Np; m++) {
                const double dr = Dr[m];
                const double ds = Ds[m];
                const double dx = dr*drdx + ds*dsdx;
                const double dy = dr*drdy + ds*dsdy;
                
                for (int fld=0; fld<Nfield; fld++) {
                    rhs[fld] += -(dx*E[sk] + dy*G[sk]);
                    sk++;
                }
            }
            // reset rhs values equal to the volume term
            for(int fld=0;fld<Nfield;fld++){
                f_rhsQ[fld][k][n] = rhs[fld];
            }
        }
    }
    
    return;
}

static void NdgSurfTerm2d
(NdgUMeshUnion *mesh,
 void (*evaluateRHSNodalFlux)(int ind, double *val, double *E, double *G, double *H),
 void (*evaluateRHSNumericalFlux)(int ind, double *fM, double *fP, double nx, double ny, double nz, double *flux),
 void (*evaluateRHSNeighbourValue)(double *fM, double *fP, double *fext, double nx, double ny, double nz, NdgEdgeType edgeNodeType, double *fextP),
 int Nfield,
 double ***f_valQ,
 double ***f_extQ,
 double ***f_rhsQ)
{

    NdgCell *cell = mesh->cell;
    const int K = mesh->K;
    const int Np = mesh->cell->Np;
    const int TNfp = mesh->cell->TNfp;
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k<K; k++) {
        
        double fluxS[Nfield]; // numerical flux
        double fM[Nfield];
        double fP[Nfield];
        double fext[Nfield];
        double fextP[Nfield];
        double dflux[Nfield*TNfp];
        
        register int sk = 0;
        for (int m=0; m<TNfp; m++) {
            const int idM = mesh->eidM[k][m] - 1;
            const int idP = mesh->eidP[k][m] - 1;
            const NdgEdgeType type = mesh->eidType[k][m];
            
            double nx = mesh->nx[k][m];
            double ny = mesh->ny[k][m];
            double Js = mesh->Js[k][m];
            
            for (int fld=0; fld<Nfield; fld++) {
                fM[fld] = f_valQ[fld][0][idM];
                fP[fld] = f_valQ[fld][0][idP];
                fext[fld] = f_extQ[fld][0][idM];
            }
            // get the external value
            evaluateRHSNeighbourValue(fM, fP, fext, nx, ny, 0, type, fextP);
            // get the numerical flux
            evaluateRHSNumericalFlux(idM, fM, fP, nx, ny, 0, fluxS);
            
            double E[Nfield], G[Nfield];            
            evaluateRHSNodalFlux(idM, fM, E, G, NULL);
            
            for (int fld=0; fld<Nfield; fld++) {
                dflux[sk++] = Js*( nx*E[fld] + ny*G[fld] - fluxS[fld] );
            }
        }
        
        // times the LIFT matrix
        for(int n=0; n<Np; ++n){
            
            const double *ptLIFT = cell->LIFTt[n];
            const double J = mesh->J[k][n];
            double rhsQ[Nfield];

            for (int fld=0; fld<Nfield; fld++) {
                rhsQ[fld] = 0.;
            }
            
            sk = 0;
            for(int m=0;m<TNfp;++m){
                const float L = ptLIFT[m];
                for (int fld=0; fld<Nfield; fld++) {
                    rhsQ[fld] += L*dflux[sk++];
                }
            }
            
            for (int fld=0; fld<Nfield; fld++) {
                f_rhsQ[fld][k][n] += rhsQ[fld]/J;
            }
        }
    }
    return;
}

/** calculate the rhs term */
void NdgRHS2d(NdgSolver *solver, NdgPhysUnion *phys, double ***f_rhsQ){
    
    NdgUMeshUnion *mesh = phys->mesh;
    int Nfield = phys->Nfield;
    double ***f_Q = phys->fvar;
    
#ifdef PROFILE
    clock_t tstart = clock()/CLOCKS_PER_SEC;
#endif
    
    NdgVolumeTerm2d(mesh,
                    solver->evaluateRHSNodalFlux,
                    Nfield, f_Q,
                    f_rhsQ);
#ifdef PROFILE
    clock_t tend = clock()/CLOCKS_PER_SEC;
    elapsedTimeVolumeTerm += tend - tstart;
    tstart = clock()/CLOCKS_PER_SEC;
#endif
    
    NdgSurfTerm2d(mesh,
                  solver->evaluateRHSNodalFlux,
                  solver->evaluateRHSNumericalFlux,
                  solver->evaluateRHSNeighbourValue,
                  Nfield, f_Q, phys->fext,
                  f_rhsQ);
#ifdef PROFILE
    tend = clock()/CLOCKS_PER_SEC;
    elapsedTimeSurfaceTerm += tend - tstart;
#endif
    
    return;
}