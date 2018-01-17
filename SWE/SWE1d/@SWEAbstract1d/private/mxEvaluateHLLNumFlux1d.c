#include <math.h>
#include "mex.h"
#include "mxSWE1d.h"

void evaluateHLLFormula(double hmin,
                        double gra,
                        double hM,
                        double qnM,
                        double hP,
                        double qnP,
                        double* Fhn,
                        double* Fqxn) {
  //   double FhM, FqnM, FhP, FqnP;
  double Em[2], Ep[2];
  evaluateSurfFluxTerm(hmin, gra, hM, qnM, Em);
  evaluateSurfFluxTerm(hmin, gra, hP, qnP, Ep);

  /* calculation of wave speed */
  double sM, sP;
  if ((hM > hmin) & (hP > hmin)) {
    double unM = qnM / hM;
    double unP = qnP / hP;
    double us = (0.5 * (unM + unP) + sqrt(gra * hM) - sqrt(gra * hP));
    double cs = (0.5 * (sqrt(gra * hM) + sqrt(gra * hP)) + 0.25 * (unM - unP));

    sM = min(unM - sqrt(gra * hM), us - cs);
    sP = max(unP + sqrt(gra * hP), us + cs);
  } else if ((hM > hmin) & (hP <= hmin)) {
    double unM = qnM / hM;
    sM = (unM - sqrt(gra * hM));
    sP = (unM + 2 * sqrt(gra * hM));
  } else if ((hM <= hmin) & (hP > hmin)) {
    double unP = qnP / hP;
    sM = (unP - 2 * sqrt(gra * hP));
    sP = (unP + sqrt(gra * hP));
  } else { /* both dry element */
    sM = 0;
    sP = 0;
  }

  /* HLL function */
  if ((sM >= 0) & (sP > 0)) {
    *Fhn = Em[0];
    *Fqxn = Em[1];
  } else if ((sM < 0) & (sP > 0)) {
    *Fhn = (sP * Em[0] - sM * Ep[0] + sM * sP * (hP - hM)) / (sP - sM);
    *Fqxn = (sP * Em[1] - sM * Ep[1] + sM * sP * (qnP - qnM)) / (sP - sM);
  } else if ((sM < 0) & (sP <= 0)) {
    *Fhn = Ep[0];
    *Fqxn = Ep[1];
  } else if ((sM == 0) & (sP == 0)) {
    *Fhn = 0;
    *Fqxn = 0;
  } else {
    mexPrintf("Matlab:%s:ErrWaveSpeed\n", __FILE__);
    mexPrintf("The wave speed computation occurs an error.");
  }
  return;
}

#define NRHS 5
#define NLHS 1

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  /* check input & output */
  if (nrhs != NRHS) {
    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NRHS);
  }

  if (nlhs != NLHS) {
    mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NLHS);
  }

  double hcrit = mxGetScalar(prhs[0]);
  double gra = mxGetScalar(prhs[1]);
  double* nx = mxGetPr(prhs[2]);
  double* fm = mxGetPr(prhs[3]);
  double* fp = mxGetPr(prhs[4]);

  const mwSize* dims = mxGetDimensions(prhs[4]);
  const size_t TNfp = dims[0];
  const size_t K = dims[1];

  const size_t NdimOut = 3;
  const mwSize dimOut[3] = {TNfp, K, NVAR};
  plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double* Fh = mxGetPr(plhs[0]);
  double* Fqx = Fh + TNfp * K;

  double* hm = fm;
  double* hp = fp;
  double* hum = fm + K * TNfp;
  double* hup = fp + K * TNfp;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif

  for (int k = 0; k < K; k++) {
    for (int n = 0; n < TNfp; n++) {
      const size_t sk = k * TNfp + n;

      double fm[2] = {hm[sk], hum[sk]};
      double fp[2] = {hp[sk], hup[sk]};

      const double nx_ = nx[sk];
      double qnM = hum[sk] * nx_;
      double qnP = hup[sk] * nx_;

      double Fhn, Fqn;
      evaluateHLLFormula(hcrit, gra, fm[0], qnM, fp[0], qnP, &Fhn, &Fqn);
      Fh[sk] = Fhn;
      Fqx[sk] = Fqn * nx_;
    }
  }
  return;
}