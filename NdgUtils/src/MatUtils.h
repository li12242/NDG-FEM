//
//  MatUtils.h
//  NDGOM_Mex
//
//  Created by li12242 on 17/10/9.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#ifndef __NDGOM_Mex__MatUtils__
#define __NDGOM_Mex__MatUtils__

#include <stdio.h>
#include <stdlib.h>

#ifdef MATLAB_MEX_FILE
//#define calloc mxCalloc
#define printf mexPrintf
#endif

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

/** allocate a vector of matrix
 @param[in] Nmat number of matrix
 @param[in] Nrow 1st dimension of the matrix
 @param[in] Ncol 2nd dimension of the matrix
 @return pointer to the matrix vector
 */
double ***mallocMatrixVector(int Nmat, int Nrow, int Ncol);

/** allocate memory for double matrix */
double **mallocMatrix(int Nrows, int Ncols);

/** allocate memory for integer matrix */
int **mallocIntMatrix(int Nrows, int Ncols);

/** free a vector of matrix */
void freeMatrixVector(double ***M);

/** free the memory of double matrix */
void freeMatrix(double **M);

/** free the memory of integer matrix */
void freeIntMatrix(int **M);

/** transpose matrix Mat of Nrows and Ncols */
double** transMatrix(int Nrows, int Ncols, double **Mat);

#endif /* defined(__NDGOM_Mex__MatUtils__) */
