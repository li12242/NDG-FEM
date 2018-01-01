#include <math.h>
#include "mex.h"
#include "mxSWE2d.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * @brief Calculation of the HLL numerical flux
 * @param [in] hmin - threadhold of the water depth
 * @param [in] gra - gravity accelerated
 * @param [in] hM, qnM, qvM - water depth on local element node
 * @param [in] hP, qnP, qvP - water depth on adjacent element node
 * @param [out] Fhn, Fqxn, Fqyn - HLL numerical flux;
 */
void evaluateHLLFormula(double hmin,
                        double gra,
                        double hM,
                        double qnM,
                        double qvM,
                        double hP,
                        double qnP,
                        double qvP,
                        double* Fhn,
                        double* Fqxn,
                        double* Fqyn) {
  double Em[3], Ep[3];
  double Gm[3], Gp[3];
  evaluateSurfFluxTerm(hmin, gra, hM, qnM, qvM, Em, Gm);
  evaluateSurfFluxTerm(hmin, gra, hP, qnP, qvP, Ep, Gp);

#if 0
  mexPrintf("hm=%f, qnM=%f, qvM=%f, Em=[%f, %f, %f], Gm=[%f, %f, %f]\n", hM,
            qnM, qvM, Em[0], Em[1], Em[2], Gm[0], Gm[1], Gm[2]);
  mexPrintf("hp=%f, qnP=%f, qvP=%f, Ep=[%f, %f, %f], Gp=[%f, %f, %f]\n", hP,
            qnP, qvP, Em[0], Em[1], Em[2], Gm[0], Gm[1], Gm[2]);
#endif
  double sM, sP, us, cs, unM, unP;
  if ((hM > hmin) & (hP > hmin)) {
    unM = qnM / hM;
    unP = qnP / hP;
    us = (double)(0.5 * (unM + unP) + sqrt(gra * hM) - sqrt(gra * hP));
    cs = (double)(0.5 * (sqrt(gra * hM) + sqrt(gra * hP)) + 0.25 * (unM - unP));

    sM = (double)min(unM - sqrt(gra * hM), us - cs);
    sP = (double)max(unP + sqrt(gra * hP), us + cs);
  } else if ((hM > hmin) & (hP <= hmin)) {
    unM = qnM / hM;
    sM = (double)(unM - sqrt(gra * hM));
    sP = (double)(unM + 2 * sqrt(gra * hM));
  } else if ((hM <= hmin) & (hP > hmin)) {
    unP = qnP / hP;
    sM = (double)(unP - 2 * sqrt(gra * hP));
    sP = (double)(unP + sqrt(gra * hP));
  } else { /* both dry element */
    sM = 0;
    sP = 0;
  }
  /* HLL flux function */
  if ((sM >= 0) & (sP > 0)) {
    *Fhn = Em[0];
    *Fqxn = Em[1];
    *Fqyn = Em[2];
  } else if ((sM < 0) & (sP > 0)) {
    *Fhn = (sP * Em[0] - sM * Ep[0] + sM * sP * (hP - hM)) / (sP - sM);
    *Fqxn = (sP * Em[1] - sM * Ep[1] + sM * sP * (qnP - qnM)) / (sP - sM);
    *Fqyn = (sP * Em[2] - sM * Ep[2] + sM * sP * (qvP - qvM)) / (sP - sM);
  } else if ((sM < 0) & (sP <= 0)) {
    *Fhn = Ep[0];
    *Fqxn = Ep[1];
    *Fqyn = Ep[2];
  } else if ((fabs(sM) < EPS) & (fabs(sP) < EPS)) {
    *Fhn = Em[0];
    *Fqxn = Em[1];
    *Fqyn = Em[2];
  } else {
    mexPrintf("Matlab:%s:ErrWaveSpeed\n", __FILE__);
    mexPrintf("The wave speed computation occurs an error.");
  }
  return;
}

#define NRHS 6
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
  double* ny = mxGetPr(prhs[3]);
  double* fm = mxGetPr(prhs[4]);
  double* fp = mxGetPr(prhs[5]);

  const mwSize* dims = mxGetDimensions(prhs[4]);
  const size_t TNfp = dims[0];
  const size_t K = dims[1];

  const size_t NdimOut = 3;
  const mwSize dimOut[3] = {TNfp, K, NVAR};
  plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double* Fh = mxGetPr(plhs[0]);
  double* Fqx = Fh + TNfp * K;
  double* Fqy = Fh + 2 * TNfp * K;

  double* hm = fm;
  double* hum = fm + K * TNfp;
  double* hvm = fm + 2 * K * TNfp;

  double* hp = fp;
  double* hup = fp + K * TNfp;
  double* hvp = fp + 2 * K * TNfp;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif

  for (int k = 0; k < K; k++) {
    for (int n = 0; n < TNfp; n++) {
      const size_t sk = k * TNfp + n;

      double fm[3] = {hm[sk], hum[sk], hvm[sk]};
      double fp[3] = {hp[sk], hup[sk], hvp[sk]};

      const double nx_ = nx[sk];
      const double ny_ = ny[sk];

      const double qnM = fm[1] * nx_ + fm[2] * ny_;
      const double qvM = -fm[1] * ny_ + fm[2] * nx_;
      const double qnP = fp[1] * nx_ + fp[2] * ny_;
      const double qvP = -fp[1] * ny_ + fp[2] * nx_;

      double Fhns, Fqns, Fqvs;
      evaluateHLLFormula(hcrit, gra, fm[0], qnM, qvM, fp[0], qnP, qvP, &Fhns,
                         &Fqns, &Fqvs);

      Fh[sk] = Fhns;
      Fqx[sk] = (Fqns * nx_ - Fqvs * ny_);
      Fqy[sk] = (Fqns * ny_ + Fqvs * nx_);

#if 0
      mexPrintf("k=%d, n=%d, hm=%f, hp=%f, hum=%f, hup=%f, hvm=%f, hvp=%f\n", k,
                n, fm[0], fp[0], qnM, qnP, qvM, qvP);
      mexPrintf("k=%d, n=%d, Fh=%f, Fqxn=%f, Fqyn=%f\n", k, n, Fhns, Fqns,
                Fqvs);
      mexPrintf("k=%d, n=%d, Fh=%f, Fqx=%f, Fqy=%f\n", k, n, Fh[sk], Fqx[sk],
                Fqy[sk]);
#endif
    }
  }

  return;
}