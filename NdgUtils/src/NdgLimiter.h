//
//  NdgLimiter.h
//  NDG-FEM
//
//  Created by li12242 on 17/10/17.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#ifndef __NDG_FEM__NdgLimiter__
#define __NDG_FEM__NdgLimiter__

#include "Utils.h"
#include "NdgPhys.h"

/** vertex-based limiter form Li (Computers and Fluids, 2017) */
void vert_limiter(NdgPhysUnion *phys, double ***val);

#endif /* defined(__NDG_FEM__NdgLimiter__) */
