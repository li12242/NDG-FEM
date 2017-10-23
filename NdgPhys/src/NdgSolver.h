//
//  NdgPhysSolver.h
//  NDG-FEM
//
//  Created by li12242 on 17/10/17.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#ifndef __NDG_FEM__NdgPhysSolver__
#define __NDG_FEM__NdgPhysSolver__

#include "NdgPhys.h"
#include "NdgPhysIO.h"

typedef struct NdgSolver{
    
    NdgPhysUnion *phys;
    NdgOutputFileInfo *outFileInfo;
    
    /**
     @brief Function pointer to temporal discrete function
     @param[in,out] phys Pointer to NdgPhysUnion class
     */
    void (*evaluateTemporalDiscrete)(struct NdgSolver *solver);
    
    /**
     @brief Function pointer for Updating the delta time through the physical field
     @param[in] phys Pointer to NdgPhysUnion class
     @param[out] dt The time interval
     */
    void (*updatePhysTimeInterval)(struct NdgPhysUnion *phys, double *dt);
    
    /**
     @brief Function pointer for updating the external values in the NdgPhysUnion class
     @param[in,out] phys NdgPhysUnion class
     @param[in] time The computation time
     @param[in] tstep The computation step
     */
    void (*updatePhysExternalInfo)(struct NdgPhysUnion *phys, double time, int tstep);
    
    /**
     @brief Function pointer to output the result
     @param[in] phys NdgPhysUnion class
     @param[in] time The computation time
     @param[in] tstep The computation step
     */
    void (*updateOutputResult)(struct NdgPhysUnion *phys, NdgOutputFileInfo *outputFileInfo, double time, int tstep);
    
    /**
     @brief Function pointer to evaluate the RHS term
     @param[in] phys NdgPhysUnion class
     @param[out] rhs The RHS term variable
     */
    void (*evaluateRHSTerm)(struct NdgSolver *solver, struct NdgPhysUnion *phys, double ***rhs);
    
    /**
     @brief Function pointer to evaluate the flux term
     @details
     This function evaluates the flux term [E, G, H] on each node.
     @param[in] ind Index of the node
     @param[in] val Pointer to the nodal values (Nfield variables)
     @param[out] E flux term on the node
     @param[out] G flux term on the node
     @param[out] G flux term on the node
     */
    void (*evaluateRHSNodalFlux)(int ind, double *val,
                                 double *E, double *G, double *H);
    
    /**
     @brief Function pointer to evaluate the numerical flux term
     @param[in] ind Index of the node
     @param[in] fM Pointer to the values on local node
     @param[in] fP Pointer to the values on adjacent node
     @param[in] nx Outward normal vector
     @param[in] ny Outward normal vector
     @param[out] flux Pointer to the numerical flux
     */
    void (*evaluateRHSNumericalFlux)(int ind,
                                     double *fM, double *fP,
                                     double nx, double ny, double nz,
                                     double *flux);
    
    /**
     @brief Function pointer to evaluate the neighbour node value
     @param[in] fM Pointer to the values on local node
     @param[in] fP Pointer to the adjacent nodel value
     @param[in] fext Pointer to the exact value
     @param[in] nx Outward normal vector
     @param[in] ny Outward normal vector
     @param[in] edgetype The edge type on node
     @param[out] fextP The exact neighbour node value
     */
    void (*evaluateRHSNeighbourValue)(double *fM, double *fP,
                                      double *fext,
                                      double nx, double ny, double nz,
                                      NdgEdgeType edgetype, double *fextP);
    
    /**
     @brief Function pointer to evaluate the source term
     @details
     The function is used to evaluates the source term after calculating
     the RHS term, and add the source values to the rhs variables.
     @param[in] solver Pointer to the physical solver
     @param[in] var Pointer to the field variable
     @param[out] rhs Pointer to the rhs term
     */
    void (*evaluateRHSSourceTerm)(struct NdgPhysUnion *phys, double ***var, double ***rhs);
    
    /**
     @brief Function pointer to evaluate the slope limiter
     @param[in] phys NdgPhysUnion class
     @param[in,out] val Pointer to the field variable
     */
    void (*evaluateLimiter)(struct NdgPhysUnion *phys, double ***val);
    
    /**
     @brief Function to evaluate the post process
     @param[in] phys NdgPhysUnion class
     @param[in,out] val Pointer to the field variable
     */
    void (*evaluatePostFunc)(struct NdgPhysUnion *phys, double ***val);
    
}NdgSolver;


/**
 @brief Setup the function pointers to solve the physical problem
 */
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
 );

void freeNdgSolver(NdgSolver *solver);

/** calculate the RHS term and store the results in NdgPhysUnion class */
void NdgRHS2d(NdgSolver *solver, NdgPhysUnion *phys, double ***rhs);

#endif /* defined(__NDG_FEM__NdgPhysSolver__) */
