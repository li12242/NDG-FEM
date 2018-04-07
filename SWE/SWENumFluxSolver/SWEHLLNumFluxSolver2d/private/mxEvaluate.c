#include "../../SWEAbstractNumFluxSolver2d/private/SWENumFlux2d.h"

/** Evaluate HLL numerical flux function */
void evaluateHLLFunc(const double hmin, ///< water threshold
                     const double gra,  ///< gravity acceleration
                     const double hM,   ///< local water depth
                     const double qnM,  ///< local flux variable
                     const double qvM,  ///< local flux variable
                     const double hP,   ///< adjacent water depth
                     const double qnP,  ///< adjacent flux
                     const double qvP,  ///< adjacent flux
                     double *Fhn,       ///< numerical flux term
                     double *Fqxn,      ///< numerical flux term
                     double *Fqyn       ///< numerical flux term
) {
  double Em[3], Ep[3];
  double Gm[3], Gp[3];
  evaluateFluxTerm2d(hmin, gra, hM, qnM, qvM, Em, Gm);
  evaluateFluxTerm2d(hmin, gra, hP, qnP, qvP, Ep, Gp);

#if 0
  mexPrintf("hm=%f, qnM=%f, qvM=%f, Em=[%f, %f, %f], Gm=[%f, %f, %f]\n", hM,
            qnM, qvM, Em[0], Em[1], Em[2], Gm[0], Gm[1], Gm[2]);
  mexPrintf("hp=%f, qnP=%f, qvP=%f, Ep=[%f, %f, %f], Gp=[%f, %f, %f]\n", hP,
            qnP, qvP, Em[0], Em[1], Em[2], Gm[0], Gm[1], Gm[2]);
#endif

  double sM, sP, unM, unP;
  if ((hM > hmin) & (hP > hmin)) {
    unM = qnM / hM;
    unP = qnP / hP;
    double us = (double)(0.5 * (unM + unP) + sqrt(gra * hM) - sqrt(gra * hP));
    double cs =
        (double)(0.5 * (sqrt(gra * hM) + sqrt(gra * hP)) + 0.25 * (unM - unP));

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
  } else if ((fabs(sM) < TOLERR) & (fabs(sP) < TOLERR)) {
    *Fhn = Em[0];
    *Fqxn = Em[1];
    *Fqyn = Em[2];
  } else {
    mexPrintf("Matlab:%s:ErrWaveSpeed\n", __FILE__);
    mexPrintf("The wave speed computation occurs an error.");
    exit(0);
  }
  return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  FluxSolver solver = ConvertInputMexVariable2d(nlhs, nrhs, plhs, prhs);

  const size_t NdimOut = 3;
  const size_t TNfp = solver.TNfp;
  const size_t K = solver.K;
  const mwSize dimOut[3] = {TNfp, K, 3};
  plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double *Fh = mxGetPr(plhs[0]);
  double *Fqx = Fh + TNfp * K;
  double *Fqy = Fh + 2 * TNfp * K;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
  for (int k = 0; k < K; k++) {
    for (int n = 0; n < TNfp; n++) {
      const size_t sk = k * TNfp + n;

      const double nx = solver.nx[sk];
      const double ny = solver.ny[sk];

      double qnM, qnP, qvM, qvP;
      RotateFluxToNormal2d(solver.huM[sk], solver.hvM[sk], nx, ny, &qnM, &qvM);
      RotateFluxToNormal2d(solver.huP[sk], solver.hvP[sk], nx, ny, &qnP, &qvP);
      double Fqns, Fqvs;
      evaluateHLLFunc(solver.hmin, solver.gra, solver.hM[sk], qnM, qvM,
                      solver.hP[sk], qnP, qvP, Fh + sk, &Fqns, &Fqvs);

      RotateNormalFluxToCoordinate2d(Fqns, Fqvs, nx, ny, Fqx + sk, Fqy + sk);
    }
  }
  return;
}