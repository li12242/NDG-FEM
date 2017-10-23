//
//  mxNdgPhys.h
//  NDG-FEM
//
//  Created by li12242 on 17/10/21.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#ifndef __NDG_FEM__mxNdgPhys__
#define __NDG_FEM__mxNdgPhys__

#include "NdgPhys.h"

/**
 @brief Create the NdgPhysUnion structure from Matlab variables
 @param[in] mxPhys Matlab variables of NdgPhys class
 @return phys Pointer to the NdgPhysUnion class
 */
NdgPhysUnion *mxGetNdgPhysUnion(const mxArray *mxPhys);


/**
 @brief Print the NdgPhysUnion structure to the Matlab command window
 @param[in] phys Pointer to the NdgPhysUnion class
 @param[in] message
 */
void mxPrintNdgPhys(NdgPhysUnion *phys, char *message);

/** 
 @brief Copy the initial field from the Matlab variable
 @param[in]
 */
void
mxCopyInitialPhysUnionField
(const mxArray *mxPhys,
 const char *fieldName,
 NdgPhysUnion *phys);

#endif /* defined(__NDG_FEM__mxNdgPhys__) */
