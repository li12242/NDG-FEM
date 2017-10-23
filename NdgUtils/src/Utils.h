#ifndef UTILS_H
#define UTILS_H

#include "mex.h"
#include "MatUtils.h"

/** read double scalar from Matlab mxArray */
double mxGetPropertyDouble(const mxArray *mxclass, char *fieldname);

/** read double vector from Matlab mxArray */
double *mxGetPropertyDoubleVector(const mxArray *mxclass, int N, char *fieldname);

/** read double matrix from Matlab mxArray */
double **mxGetPropertyDoubleMatrix(const mxArray *mxCell, int Nrow, int Ncol, char *fieldname);

/** read int8 type filed value from Matlab mxArray */
signed char mxGetPropertyInt8(const mxArray *mxclass, char *fieldname);

/** read integer vector from Matlab mxArray */
int *mxGetPropertyIntVector(const mxArray *mxclass, int N, char *fieldname);

/** read integer matrix from Matlab mxArray */
int **mxGetPropertyIntMatrix(const mxArray *mxEdge, int Nrow, int Ncol, char *fldname);

/** print integer matrix to Matlab command window */
void mexPrintIntMatrix(char *str, int **Mat, int Nrow, int Ncol);

/** print double matrix to Matlab command window */
void mexPrintDoubleMatrix(char *str, double **Mat, int Nrow, int Ncol);

#endif