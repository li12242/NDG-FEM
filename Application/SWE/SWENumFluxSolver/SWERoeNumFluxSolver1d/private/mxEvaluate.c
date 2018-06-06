#include "../../SWEAbstractNumFluxSolver1d/private/SWENumFlux1d.h"

// #define DEBUG

typedef struct {
  double h;
  double u;
  double c;
} RoeState;

inline void evaluateVelocity(const double hcrit, ///< depth threshold
                             const double h,     ///< depth
                             const double hu,    ///< water flux
                             double *u           ///< result velocity
) {
  if (h > hcrit) {
    *u = hu / h;
  } else {
    *u = 0.0;
  }
  return;
}

inline void evaluateRoeAverage(const double gra,   ///< gravity acceleration
                               const double hcrit, ///< water depth threshold
                               const double hM,    ///< local water depth
                               const double huM,   ///< local flux
                               const double hP,    ///< neighbour water depth
                               const double huP,   ///< neighbour flux
                               RoeState *roe       ///< averaged Roe state
) {
  double hsqrtM = sqrt(hM);
  double hsqrtP = sqrt(hP);
  double uM, uP;
  evaluateVelocity(hcrit, hM, huM, &uM);
  evaluateVelocity(hcrit, hP, huP, &uP);
  roe->h = hsqrtM * hsqrtP;
  roe->u = (uM * hsqrtM + uP * hsqrtP) / (hsqrtM + hsqrtP);
  roe->c = sqrt(gra * (hM + hP) * 0.5);
#ifdef DEBUG
  mexPrintf("Roe averaged states\n");
  mexPrintf("h = %f\nu = %f\nc = %f\n", roe->h, roe->u, roe->c);
#endif
  return;
}

inline void
evaluateRoeWaveStrength(const double hcrit,  ///< water depth threshold
                        const double hM,     ///< local water depth
                        const double huM,    ///< local flux
                        const double hP,     ///< neighbour water depth
                        const double huP,    ///< neighbour flux
                        const double nx,     ///< outward normal vector
                        const RoeState *roe, ///< roe averaged states
                        double *alpha        ///< wave strength
) {
  const double qnM = huM * nx;
  const double qnP = huP * nx;
  double unM, unP;
  evaluateVelocity(hcrit, hM, qnM, &unM);
  evaluateVelocity(hcrit, hP, qnP, &unP);
  alpha[0] = 0.5 * (hP - hM - roe->h / roe->c * (unP - unM));
  alpha[1] = 0.5 * (hP - hM + roe->h / roe->c * (unP - unM));
#ifdef DEBUG
  mexPrintf("Wave strength\n");
  mexPrintf("local velocity %f\n", unM);
  mexPrintf("neigh velocity %f\n", unP);
  mexPrintf("alpha = [%f, %f]\n", alpha[0], alpha[1]);
#endif
  return;
}

void evaluateRoeSolver(const double hmin,   ///< water depth threshold
                       const double gra,    ///< gravity acceleration
                       const double hM,     ///< local water depth
                       const double huM,    ///< local flux
                       const double hP,     ///< neighbour water depth
                       const double huP,    ///< neighbour flux
                       const double nx,     ///< outward normal vector
                       const RoeState *roe, ///< roe averaged states
                       double *Fh,          ///< roe flux on h
                       double *Fhu          ///< roe flux on hu
) {
  double EM[2];
  evaluateFluxTerm1d(hmin, gra, hM, huM, EM);
  Fh[0] = EM[0] * nx;
  Fhu[0] = EM[1] * nx;
  evaluateFluxTerm1d(hmin, gra, hP, huP, EM);
  Fh[0] += EM[0] * nx;
  Fhu[0] += EM[1] * nx;

  double alpha[2];
  evaluateRoeWaveStrength(hmin, hM, huM, hP, huP, nx, roe, alpha);
  const double unroe = roe->u * nx;
  const double lambda1 = fabs(unroe - roe->c);
  const double lambda3 = fabs(unroe + roe->c);

#ifdef DEBUG
  mexPrintf("eigenvalue lambda = [%f, %f]\n", lambda1, lambda3);
#endif
  Fh[0] -= lambda1 * alpha[0];
  Fhu[0] -= lambda1 * alpha[0] * (roe->u - roe->c * nx);

  Fh[0] -= lambda3 * alpha[1];
  Fhu[0] -= lambda3 * alpha[1] * (roe->u + roe->c * nx);

  Fh[0] *= 0.5;
  Fhu[0] *= 0.5;

#ifdef DEBUG
  mexPrintf("Roe flux = [%f, %f]\n", Fh[0], Fhu[0]);
#endif
  return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  FluxSolver solver = ConvertInputMexVariable1d(nlhs, nrhs, prhs);

  const size_t NdimOut = 3;
  const size_t K = solver.K;
  const size_t TNfp = solver.TNfp;
  const mwSize dimOut[3] = {TNfp, K, 2};
  plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double *Fh = mxGetPr(plhs[0]);
  double *Fqx = Fh + solver.TNfp * solver.K;

#ifndef DEBUG
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
#endif
  for (int k = 0; k < K; k++) {
    for (int n = 0; n < TNfp; n++) {
      const size_t sk = k * TNfp + n;
      const double hM = solver.hM[sk];
      const double hP = solver.hP[sk];
      const double huM = solver.huM[sk];
      const double huP = solver.huP[sk];
#ifdef DEBUG
      mexPrintf("k = %d, sk = %d\n", k, n);
      mexPrintf("h = [%f, %f]\nhu = [%f, %f]\n", hM, hP, huM, huP);
#endif
      if ((hM > solver.hmin) || (hP > solver.hmin)) {
        const double nx = solver.nx[sk];
        RoeState roe;
        evaluateRoeAverage(solver.gra, solver.hmin, hM, huM, hP, huP, &roe);
        evaluateRoeSolver(solver.hmin, solver.gra, hM, huM, hP, huP, nx, &roe,
                          Fh + sk, Fqx + sk);
      }
    }
  }
  return;
}