//
//  MatUtils.c
//  NDGOM_Mex
//
//  Created by li12242 on 17/10/9.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#include "MatUtils.h"

/** allocate a vector of double matrix */
double ***mallocMatrixVector(int Nmat, int Nrow, int Ncol){
    
    double ***A = (double ***)calloc(Nmat, sizeof(double **));
    A[0] = (double **)calloc(Nmat*Nrow, sizeof(double *));
    A[0][0] = (double *)calloc(Nmat*Nrow*Ncol, sizeof(double));
    for (int m=1; m<Nmat; m++) {
        A[m] = A[m-1] + Nrow;
    }
    
    for (int m = 0; m<Nmat; m++) {
        for (int n=1; n<Nrow; n++) {
            A[m][n] = A[m][n-1] + Ncol;
        }
    }
    
    return A;
}

/** transpose matrix Mat of Nrows and Ncols */
double** transMatrix(int Nrows, int Ncols, double **Mat){
    double **A = mallocMatrix(Ncols, Nrows);
    for(int n = 0; n<Ncols; n++)
        for(int m =0; m<Nrows; m++)
            A[n][m] = Mat[m][n];
    return A;
}

/** allocate column major storage for a 2D double matrix. */
double **mallocMatrix(int Nrows, int Ncols)
{
    double **A = (double **)calloc(Nrows, sizeof(double *));
    A[0] = (double *)calloc(Nrows * Ncols, sizeof(double));
    for (int n = 1; n < Nrows; ++n)
        A[n] = A[n - 1] + Ncols;
    return A;
}

/** allocate column major storage for a 2D integer matrix. */
int **mallocIntMatrix(int Nrows, int Ncols)
{
    int **A = (int **)calloc(Nrows, sizeof(int *));
    A[0] = (int *)calloc(Nrows * Ncols, sizeof(int));
    for (int n = 1; n < Nrows; ++n)
        A[n] = A[n - 1] + Ncols;
    return A;
}

/** free a vector of double matrix */
void freeMatrixVector(double ***M){
    free(M[0][0]);
    free(M[0]);
    free(M);
}

/** free 2D double matrix. */
void freeMatrix(double **M)
{
    free(M[0]);
    free(M);
    return;
}

/** free 2D integer matrix */
void freeIntMatrix(int **M)
{
    free(M[0]);
    free(M);
    return;
}