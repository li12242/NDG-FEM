//
//  UMeshIO.h
//  NDGOM_Mex
//
//  Created by li12242 on 17/10/10.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#ifndef __NDGOM_Mex__UMeshIO__
#define __NDGOM_Mex__UMeshIO__

#include "NdgMesh.h"
#include "NdgNetCDF.h"

/** print the NdgMeshUnion structure to Matlab command window. */
void mxPrintUMeshUnion(NdgUMeshUnion *mesh, char *message);

/** copy the NdgMeshUnion structure from Matlab variable. */
NdgUMeshUnion *mxGetNdgUMeshUnion(const mxArray *mxMesh);

/** 
 @brief Create the output NetCDF file and stores some basic variables into the file.
 @details
 Create a NetCDF file for output the basic mesh information and return the NdgNcFile 
 pointer. The basic mesh information includes the information of cells and vertices:
 * * EToV - vertex index in each cell
 * * EToB - edge type of faces in each cell
 * * EToR - cell types of each cell
 * * vx, vy, vz - vertex coordinates
 
 The output file also include the output variables, 'time' and 'fvar'. The 'time' variable
 contains one unlimited dimension Time, while the 'fvar' has four dimensions
 [Time, Nfield, Nelement, Npoint] and iterates from right to left.
 
 @param[in] mesh pointer of the UMeshUnion object
 @param[in] casename case name of the output NetCDF file
 @param[in] Nfield Number of physical field
 @return ncfile pointer to a new NdgNcFile variable of the output file.
 */
NdgNcFile* setupBasicUMeshUnionOutputFile(NdgUMeshUnion mesh, char *casename, int Nfield);

/** */
void freeBasicUMeshUnionOutputFile(NdgNcFile *ncfile);

/**
 Free the memory of the 
 */
void freeUMeshUnionOutpufFile(NdgNcFile *ncfile);

/** 
 @brief output the result to the NetCDF file
 @details
 */
void outputBasicUMeshUnionResult(NdgNcFile ncfile, int tstep, double time, double ***fvar);

#endif /* defined(__NDGOM_Mex__UMeshIO__) */
