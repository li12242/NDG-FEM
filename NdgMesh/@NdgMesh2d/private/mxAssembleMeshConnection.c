#ifdef _OPENMP
#include <omp.h>
#endif

#include "mex.h"
#include "math.h"

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

#define NRHS 13
#define NLHS 3

/**
 * @brief
 * @details
 * [ EToM, EToE, EToF, EToB ] = mxAssembleMeshConnection( ...
    MeshId1, MeshId2, K1, K2, Nface1, Nface2, ...
    FToV1, FToV2, EToV1, EToV2, EToM, EToE, EToF, EToB )
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* check input & output */
  if (nrhs != NRHS) {
    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NRHS);
  }

  if (nlhs != NLHS) {
    mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NLHS);
  }

  double MeshId1 = mxGetScalar(prhs[0]);
  double MeshId2 = mxGetScalar(prhs[1]);
  int K1 = (int)mxGetScalar(prhs[2]);
  int K2 = (int)mxGetScalar(prhs[3]);
  int Nface1 = (int)mxGetScalar(prhs[4]);
  int Nface2 = (int)mxGetScalar(prhs[5]);
  double *FToV1 = mxGetPr(prhs[6]);
  double *FToV2 = mxGetPr(prhs[7]);
  double *EToV1 = mxGetPr(prhs[8]);
  double *EToV2 = mxGetPr(prhs[9]);
  double *EToM = mxGetPr(prhs[10]);
  double *EToE = mxGetPr(prhs[11]);
  double *EToF = mxGetPr(prhs[12]);

  const int Nv1 = mxGetM(prhs[8]);
  const int Nv2 = mxGetM(prhs[9]);

  //   mexPrintf("Nv1 = %d, Nv2 = %d\n", Nv1, Nv2);

  const size_t ndimOut = 2;
  const mwSize dimOut[2] = {Nface1, K1};
  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
  plhs[2] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
  double *EToM1 = mxGetPr(plhs[0]);
  double *EToE1 = mxGetPr(plhs[1]);
  double *EToF1 = mxGetPr(plhs[2]);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif

  for (int k1 = 0; k1 < K1; k1++) {
    for (int f1 = 0; f1 < Nface1; f1++) {
      const int sk = f1 + k1 * Nface1;
      EToM1[sk] = MeshId1;
      EToE1[sk] = EToE[sk];
      EToF1[sk] = EToF[sk];

      const int ind = (int)EToE[sk] - 1;

      if (ind == k1) {
        // mexPrintf("k1 = %d, f1 = %d \n", k1, f1);

        int v11 = EToV1[(int)FToV1[f1 * 2] - 1 + k1 * Nv1];
        int v12 = EToV1[(int)FToV1[f1 * 2 + 1] - 1 + k1 * Nv1];
        int va1 = max(v11, v12);
        int vb1 = min(v11, v12);

        for (int k2 = 0; k2 < K2; k2++) {
          for (int f2 = 0; f2 < Nface2; f2++) {
            int v21 = EToV2[(int)FToV2[f2 * 2] - 1 + k2 * Nv2];
            int v22 = EToV2[(int)FToV2[f2 * 2 + 1] - 1 + k2 * Nv2];
            int va2 = max(v21, v22);
            int vb2 = min(v21, v22);

            if ((va1 == va2) && (vb1 == vb2)) {
              EToM1[sk] = MeshId2;
              EToE1[sk] = k2 + 1;  // start from 1
              EToF1[sk] = f2 + 1;  // start from 1
            }
          }
        }
      }
    }
  }
}