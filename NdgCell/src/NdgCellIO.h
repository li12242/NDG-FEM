//
//  NdgCellIO.h
//  NDG-FEM
//
//  Created by li12242 on 17/10/21.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#ifndef __NDG_FEM__NdgCellIO__
#define __NDG_FEM__NdgCellIO__

#include "NdgCell.h"

/** copy the NdgCell form Matlab variable to C structure. */
NdgCell *mxGetNdgCell(const mxArray *cell_mxArray);

/** print the NdgCell structure to Matlab command window. */
void mxPrintNdgCellInfo(NdgCell *cell, char *message);

#endif /* defined(__NDG_FEM__NdgCellIO__) */
