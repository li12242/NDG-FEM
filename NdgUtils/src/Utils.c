#include "Utils.h"

/**
 * @brief
 */
double mxGetPropertyDouble(const mxArray *mxclass, char *fieldname)
{
    mxArray *mxField = mxGetProperty(mxclass, 0, fieldname);
    mxArray *mxDouble;
    mexCallMATLAB(1, &mxDouble, 1, &mxField, "double");
    double dval = *(mxGetPr(mxDouble));
    mxDestroyArray(mxField);
    mxDestroyArray(mxDouble);
    return dval;
}

/**
 * @brief
 */
double *mxGetPropertyDoubleVector(const mxArray *mxclass, int N, char *fieldname)
{
    mxArray *mxField = mxGetProperty(mxclass, 0, fieldname);
    double *field_ptr = mxGetPr(mxField);
    double *dvec = (double *)calloc(N, sizeof(double));
    for (int n = 0; n < N; n++)
        dvec[n] = field_ptr[n];
    mxDestroyArray(mxField);
    return dvec;
}

/**
 * @brief
 */
double **mxGetPropertyDoubleMatrix(const mxArray *mxCell, int Nrow, int Ncol, char *fieldname)
{
    double **Mat = mallocMatrix(Nrow, Ncol);
    double *dTemp = mxGetPropertyDoubleVector(mxCell, Nrow * Ncol, fieldname);
    int sk = 0;
    for (int n = 0; n < Nrow; n++)
    {
        for (int m = 0; m < Ncol; m++)
        {
            Mat[n][m] = dTemp[sk];
            sk++;
        }
    }
    free(dTemp);
    return Mat;
}

/**
 * @brief
 */
int *mxGetPropertyIntVector(const mxArray *mxclass, int N, char *fieldname)
{
    mxArray *mxField = mxGetProperty(mxclass, 0, fieldname);
    double *field_ptr = mxGetPr(mxField);
    int *ivec = (int *)calloc(N, sizeof(int));
    for (int n = 0; n < N; n++)
        ivec[n] = (int)field_ptr[n];
    
    mxDestroyArray(mxField);
    return ivec;
}

/** @brief */
int **mxGetPropertyIntMatrix(const mxArray *mxEdge, int Nrow, int Ncol, char *fldname)
{
    int **intMat = mallocIntMatrix(Nrow, Ncol);
    int *intVec = mxGetPropertyIntVector(mxEdge, Nrow * Ncol, fldname);
    int sk = 0;
    for (int n = 0; n < Nrow; n++)
    {
        for (int m = 0; m < Ncol; m++)
        {
            intMat[n][m] = intVec[sk];
            sk++;
        }
    }
    free(intVec);
    return intMat;
}

/** @brief */
void mexPrintIntMatrix(char *str, int **Mat, int Nrow, int Ncol){
    mexPrintf("%s\n", str);
    for (int n = 0; n < Nrow; n++){
        mexPrintf("\t\t ");
        for (int m = 0; m < Ncol; m++)
            mexPrintf("%d, ", Mat[n][m]);
        mexPrintf("\n");
    }
}

/** @brief */
void mexPrintDoubleMatrix(char *str, double **Mat, int Nrow, int Ncol){
    mexPrintf("%s\n", str);
    for (int n = 0; n < Nrow; n++){
        mexPrintf("\t\t ");
        for (int m = 0; m < Ncol; m++)
            mexPrintf("%f, ", Mat[n][m]);
        mexPrintf("\n");
    }
}