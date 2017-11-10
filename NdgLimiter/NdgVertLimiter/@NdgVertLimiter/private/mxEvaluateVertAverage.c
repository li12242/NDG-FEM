//
//  main.c
//  matEvaluateVertAverage
//
//  Created by li12242 on 17/11/1.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#include "mex.h"
#define INF 10.0e9

#ifdef _OPENMP
#include <omp.h>
#endif

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* check input & output */
    if (nrhs != 6){
        mexErrMsgIdAndTxt("Matlab:mxEvaluateVertAverage:InvalidNumberInput",
                          "4 inputs required.");
    }
    
    if (nlhs != 3){
        mexErrMsgIdAndTxt("Matlab:mxEvaluateVertAverage:InvalidNumberOutput",
                          "1 output required.");
    }
    
    /* get inputs */
    //double *cvar = mxGetPr(prhs[0]);
    int Nv = (int)mxGetScalar(prhs[1]);
    double *NkToV = mxGetPr(prhs[2]); // number of element connecting to each vertex
    double *VToM = mxGetPr(prhs[3]);
    double *VToK = mxGetPr(prhs[4]);
    double *VToW = mxGetPr(prhs[5]);
    
    /* get dimensions */
    size_t maxNk = mxGetM(prhs[3]); // maximum number of cells connecting to one vertex
//    size_t K = mxGetN(prhs[2]);  // number of elements
    
    /* allocate output array */
    plhs[0] = mxCreateDoubleMatrix((mwSize)Nv, (mwSize)1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)Nv, (mwSize)1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix((mwSize)Nv, (mwSize)1, mxREAL);
    double *fvert = mxGetPr(plhs[0]);
    double *fvmin = mxGetPr(plhs[1]);
    double *fvmax = mxGetPr(plhs[2]);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int n=0; n<Nv; n++) {
        int Nk = NkToV[n]; // number of cells connecting to vertex n
        fvmax[n] = -INF;
        fvmin[n] = INF;
        for (int k=0; k<Nk; k++) {
            int meshId = (int) VToM[n*maxNk + k] - 1;
            int cellId = (int) VToK[n*maxNk + k] - 1;
            double w = VToW[n*maxNk + k];
            mxArray *mcvar = mxGetCell(prhs[0], meshId);
            if (mcvar == NULL) {
                mexErrMsgIdAndTxt("Matlab:mxEvaluateVertAverage:AccessToMeshField",
                                  "Access to the mesh field failed.");
            }
            double *cvar = mxGetPr(mcvar);
            double temp = cvar[cellId];
            fvert[n] += w*temp;
            fvmax[n] = max(fvmax[n], temp);
            fvmin[n] = min(fvmin[n], temp);
        }
    }
    return;
}
