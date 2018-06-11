#include "../../SWEAbstractNumFluxSolver1d/private/SWENumFlux1d.h"

void evaluateHLLFormula(const double hmin, ///< water threshold
                        const double gra,  ///< gravity acceleration
                        const double hM,   ///< local water depth
                        const double qnM,  ///< local flux variable
                        const double hP,   ///< adjacent water depth
                        const double qnP,  ///< adjacent flux
                        double *Fhn,       ///< numerical flux term
                        double *Fqxn       ///< numerical flux term
) {
  //   double FhM, FqnM, FhP, FqnP;
  double Em[2], Ep[2];
  evaluateFluxTerm1d(hmin, gra, hM, qnM, Em);
  evaluateFluxTerm1d(hmin, gra, hP, qnP, Ep);

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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  FluxSolver solver = ConvertInputMexVariable1d(nlhs, nrhs, prhs);

  const size_t NdimOut = 3;
  const mwSize dimOut[3] = {solver.TNfp, solver.K, 2};
  plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double *Fh = mxGetPr(plhs[0]);
  double *Fqx = Fh + solver.TNfp * solver.K;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
  for (int k = 0; k < solver.K; k++) {
    for (int n = 0; n < solver.TNfp; n++) {
      const size_t sk = k * solver.TNfp + n;

      const double nx = solver.nx[sk];

      double qnM, qnP;
      RotateFluxToNormal1d(solver.huM[sk], nx, &qnM);
      RotateFluxToNormal1d(solver.huP[sk], nx, &qnP);
      double Fqns;
      evaluateHLLFormula(solver.hmin, solver.gra, solver.hM[sk], qnM,
                         solver.hP[sk], qnP, Fh + sk, &Fqns);

      RotateNormalFluxToCoordinate1d(Fqns, nx, Fqx + sk);
    }
  }
  return;
}