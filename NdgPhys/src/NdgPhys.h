//
//  NdgPhys.h
//  NDG-FEM
//
//  Created by li12242 on 17/10/16.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#ifndef __NDG_FEM__NdgPhys__
#define __NDG_FEM__NdgPhys__

#include "NdgMesh.h"

/** 
 @brief Enumerations for updating the open boundary types
 
 @details
 The boundary condition can be updated by analytical function or NetCDF files.
 If no external informations is required, the type should be NdgBCNone.
 */
typedef enum {
    NdgBCNone = 0, ///< none external information
    NdgBCFunc = 1, ///< external information is obtained from analytical function
    NdgBCFile = 2, ///< external information is obtained from NetCDF file
}NdgBCType;

/**
 @brief Enumeration for the interval type
 
 @details
 The NdgIntervalType determines the way to update the external information
 from the open boundary, output the results, or the time interval.
 */
typedef enum {
    NdgIntervalConst = 1, ///< const time interval for time discrete scheme
    NdgIntervalDeltaTime = 2, ///< update or output by delta time
    NdgIntervalDeltaStep = 3, ///< update or output by delta step
}NdgIntervalType;

/** 
 @brief Enumerations for the slope limiter type
 */
typedef enum {
    NdgLimiterNone = 0, ///< none limiter
    NdgLimiterVert = 1, ///< vertex-based limiter
}NdgLimiterType;

/** 
 @brief Enumerations for the temporal discrete scheme
 */
typedef enum {
    NdgTemporalDiscreteEuler = 1, ///< Euler scheme
    NdgTemporalDiscreteRK45 = 2, ///< 4 order 5 stages SSP-RK scheme
}NdgTemporalDiscreteType;

/** 
 @brief PhysUnion structure 
 */
typedef struct NdgPhysUnion{
    
    int Nfield; ///< number of physical filed
    NdgUMeshUnion *mesh; ///< mesh object
    double finalTime; ///< final times
    
    NdgLimiterType limiterType; ///< enumeration for limiter type
    NdgTemporalDiscreteType temporalDiscreteType; ///< enumeration for temporal discrete type
    NdgIntervalType timeIntervalType; ///< enumeration for evaluating delta time
    double timeInterval; ///< delta time
    
    NdgBCType obcType; ///< enumeration for open boundary type
    NdgIntervalType obcIntervalType; ///< enumeration for time interval updating external time
    double obcTimeInterval; ///< time interval for updating the obc (if obcIntervalType == NdgIntervalDeltaTime)
    int obcStepInterval; ///< step interval for updating the obc (if obcIntervalType == NdgIntervalDeltaStep)
    char obcFileName[NC_MAX_NAME_LEN]; ///< the open boundary file name
    
    NdgIntervalType outputIntervalType;
    char outputNetcdfCaseName[NC_MAX_NAME_LEN];
    double outputTimeInterval;
    int outputStepInterval;
    
    double ***fvar;
    double ***fext;
    double ***frhs;
    
}NdgPhysUnion;

/** 
 @brief Free the NdgPhysUnion structure variable 
 @param[in] phys Free the memory of NdgPhysUnion class
 */
void freeNdgPhys(NdgPhysUnion *phys);

#include "mxNdgPhys.h"

#endif /* defined(__NDG_FEM__NdgPhys__) */
