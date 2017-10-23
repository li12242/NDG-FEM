//
//  NdgTerporalDiscrete.c
//  NDG-FEM
//
//  Created by li12242 on 17/10/17.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#include "NdgTerporalDiscrete.h"

#ifdef _OPENMP
#include "libiomp/omp.h"
#endif

void Euler(NdgSolver *solver){
    return;
}

void SSP_RK45(NdgSolver *solver)
{
    double time = 0, dt = 0;
    int tstep = 0;
    double rk4a[5], rk4b[5], rk4c[6];
    
    rk4a[0] =              0.0;
    rk4a[1] =  -567301805773.0 / 1357537059087.0;
    rk4a[2] = -2404267990393.0 / 2016746695238.0;
    rk4a[3] = -3550918686646.0 / 2091501179385.0;
    rk4a[4] = -1275806237668.0 /  842570457699.0;
    rk4b[0] =  1432997174477.0 /  9575080441755.0;
    rk4b[1] =  5161836677717.0 / 13612068292357.0;
    rk4b[2] =  1720146321549.0 /  2090206949498.0;
    rk4b[3] =  3134564353537.0 /  4481467310338.0;
    rk4b[4] =  2277821191437.0 / 14882151754819.0;
    rk4c[0] =              0.0;
    rk4c[1] =  1432997174477.0 / 9575080441755.0;
    rk4c[2] =  2526269341429.0 / 6820363962896.0;
    rk4c[3] =  2006345519317.0 / 3224310063776.0;
    rk4c[4] =  2802321613138.0 / 2924317926251.0;
    rk4c[5] =              1.0;
    
    NdgPhysUnion *phys = solver->phys;
    const int Nfield = phys->Nfield;
    const int K = phys->mesh->K;
    const int Np = phys->mesh->cell->Np;
    const double ftime = phys->finalTime;
    
    double ***resQ = mallocMatrixVector(Nfield, K, Np);
    double ***rhsQ = phys->frhs;
    double ***f_Q = phys->fvar;
    
    while (time < ftime) {
        solver->updatePhysTimeInterval(phys, &dt);
        if (time + dt > ftime) { dt = ftime - time; }
        
        for (int interk=0; interk<5; interk++) {
            const double frka = rk4a[interk];
            const double frkb = rk4b[interk];
            const double frkc = rk4c[interk];
            
            solver->updatePhysExternalInfo(phys, time + frkc*dt, tstep);
            
            solver->evaluateRHSTerm(solver, phys, rhsQ);
            solver->evaluateRHSSourceTerm(phys, f_Q, rhsQ);
            
            for (int fld=0; fld<Nfield; fld++) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
                for (int n=0; n<K*Np; n++) {
                    resQ[fld][0][n] = frka*resQ[fld][0][n] + dt*rhsQ[fld][0][n];
                    f_Q[fld][0][n] += frkb*resQ[fld][0][n];
                }
            }
            
            solver->evaluateLimiter(phys, f_Q);
            solver->evaluatePostFunc(phys, f_Q);
        }
        time += dt;
#ifdef PROFILE
        //mexPrintf("finishing %f...\n", time/ftime);
#endif
        solver->updateOutputResult(phys, solver->outFileInfo, time, tstep);
    }

    freeMatrixVector(resQ);
    return;
}