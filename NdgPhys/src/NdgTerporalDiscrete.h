//
//  NdgTerporalDiscrete.h
//  NDG-FEM
//
//  Created by li12242 on 17/10/17.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#ifndef __NDG_FEM__NdgTerporalDiscrete__
#define __NDG_FEM__NdgTerporalDiscrete__

#include "NdgSolver.h"

/** 1 step Euler advance scheme */
void Euler(NdgSolver *solver);

/** 4 order 5 stages strong stability preserving Runge-Kutta scheme */
void SSP_RK45(NdgSolver *solver);

#endif /* defined(__NDG_FEM__NdgTerporalDiscrete__) */
