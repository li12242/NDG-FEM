#ifdef _OPENMP
#include <omp.h>
#endif

#include "mex.h"
#include "blas.h"

#if !defined(_WIN32)
#define dgemm dgemm_
#endif
// #if !defined(_WIN32)
// #define dgemm dgemm_
// #endif

#define DEBUG 0

#define NRHS 10
#define NLHS 1

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

  double *invM = mxGetPr(prhs[0]);
  double *Mb = mxGetPr(prhs[1]);
  double *FToE = mxGetPr(prhs[2]);
  double *FToN1 = mxGetPr(prhs[3]);
  double *FToN2 = mxGetPr(prhs[4]);
  double *Js = mxGetPr(prhs[5]);
  double *J = mxGetPr(prhs[6]);
  double *fluxM = mxGetPr(prhs[7]);
  double *fluxP = mxGetPr(prhs[8]);
  double *fluxS = mxGetPr(prhs[9]);

  // dims = mxGetDimensions(prhs[6]);
  const int Np = mxGetM(prhs[6]);  // num of interp nodes
  const int K = mxGetN(prhs[6]);   // num of elements
  const mwSize *dims = mxGetDimensions(prhs[7]);
  const int Nfp = dims[0];
  const int Ne = dims[1];  // num of edges
  int Nfield;

  if (mxGetNumberOfDimensions(prhs[7]) > 2) {
    Nfield = dims[2];
  } else {
    Nfield = 1;  // fluxM is a 2D matrix
  }

  const size_t ndimOut = 3;
  const mwSize dimOut[3] = {Np, K, Nfield};

  plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
  double *frhs = mxGetPr(plhs[0]);

  char *chn = "N";
  double one = 1.0, zero = 0.0;
  ptrdiff_t oneI = 1;
  ptrdiff_t np = Np;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
  for (int fld = 0; fld < Nfield; fld++) {
    double *rhs = frhs + Np * K * fld;
    double *fluxM_ = fluxM + Nfp * Ne * fld;
    double *fluxP_ = fluxP + Nfp * Ne * fld;
    double *fluxS_ = fluxS + Nfp * Ne * fld;

    for (int k = 0; k < Ne; k++) {  // evaluate rhs on each edge
      const int e1 = (int)FToE[2 * k] - 1;
      const int e2 = (int)FToE[2 * k + 1] - 1;
      const int ind1 = e1 * Np - 1;
      const int ind2 = e2 * Np - 1;

      const int ind = k * Nfp;
      double rhsM[Nfp], rhsP[Nfp];
      for (int n = 0; n < Nfp; n++) {
        rhsM[n] = 0;
        rhsP[n] = 0;
      }

      for (int n = 0; n < Nfp; n++) {
        const int sk = n + ind;
        double dfM = fluxM_[sk] - fluxS_[sk];
        double dfP = fluxP_[sk] - fluxS_[sk];
        double j = Js[sk];

        double *mb = Mb + n * Nfp;

        for (int m = 0; m < Nfp; m++) {
          rhsM[m] += mb[m] * j * dfM;
          rhsP[m] -= mb[m] * j * dfP;
        }
      }

      for (int n = 0; n < Nfp; n++) {
        const int sk = n + ind;
        const int m1 = (int)FToN1[sk] + ind1;
        const int m2 = (int)FToN2[sk] + ind2;
        rhs[m1] += rhsM[n];
        rhs[m2] += rhsP[n];
      }
    }
    double temp[Np];
    for (int k = 0; k < K; k++) {
      double *rhs_ = rhs + k * Np;
      double *j = J + k * Np;

      dgemm(chn, chn, &np, &oneI, &np, &one, invM, &np, rhs_, &np, &zero, temp,
            &np);
      // copy rhs
      for (int n = 0; n < Np; n++) {
        rhs_[n] = temp[n] / j[n];
      }
    }
  }
  return;
}