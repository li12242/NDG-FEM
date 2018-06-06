#include "../../SWEAbstractNumFluxSolver2d/private/SWENumFlux2d.h"

// #define DEBUG

typedef struct {
  double h;
  double u;
  double v;
  double c;
} RoeState;

inline void evaluateVelocity(const double hcrit, ///< depth threshold
                             const double h,     ///< depth
                             const double hu,    ///< water flux
                             const double hv,    ///< water flux
                             double *u,          ///< result velocity
                             double *v           ///< velocity
) {
  if (h > hcrit) {
    *u = hu / h;
    *v = hv / h;
  } else {
    *u = 0.0;
    *v = 0.0;
  }
  return;
}

inline void evaluateRoeAverage(const double gra,   ///< gravity acceleration
                               const double hcrit, ///< water depth threshold
                               const double hM,    ///< local water depth
                               const double huM,   ///< local flux
                               const double hvM,   ///< local flux
                               const double hP,    ///< neighbour water depth
                               const double huP,   ///< neighbour flux
                               const double hvP,   ///< neighbour flux
                               RoeState *roe       ///< averaged Roe state
) {
  //   double uM, vM;
  //   double uP, vP;
  //   evaluateVelocity(hcrit, hM, huM, hvM, &uM, &vM);
  //   evaluateVelocity(hcrit, hP, huP, hvP, &uP, &vP);
  double hsqrtM = sqrt(hM);
  double hsqrtP = sqrt(hP);
  roe->h = hsqrtM * hsqrtP;
  double uM, uP;
  double vM, vP;
  evaluateVelocity(hcrit, hM, huM, hvM, &uM, &vM);
  evaluateVelocity(hcrit, hP, huP, hvP, &uP, &vP);
  roe->u = (uM * hsqrtM + uP * hsqrtP) / (hsqrtM + hsqrtP);
  roe->v = (vM * hsqrtM + vP * hsqrtP) / (hsqrtM + hsqrtP);
  roe->c = sqrt(gra * (hM + hP) * 0.5);
#ifdef DEBUG
  mexPrintf("Roe averaged states\n");
  mexPrintf("h = %f\nu = %f\nv = %f\nc = %f\n", roe->h, roe->u, roe->v, roe->c);
#endif
  return;
}

void evaluateRoeWaveStrength(const double hcrit,  ///< water depth threshold
                             const double hM,     ///< local water depth
                             const double huM,    ///< local flux
                             const double hvM,    ///< local flux
                             const double hP,     ///< neighbour water depth
                             const double huP,    ///< neighbour flux
                             const double hvP,    ///< neighbour flux
                             const double nx,     ///< outward normal vector
                             const double ny,     ///< outward normal vector
                             const RoeState *roe, ///< roe averaged states
                             double *alpha        ///< wave strength
) {
  //   const double dh = hP - hM;
  //   alpha[0] = dh * 0.5 + (roe->u * dh - (huP - huM)) * 0.5 / roe->c;
  //   alpha[1] = hvP - hvM - roe->v * dh;
  //   alpha[2] = dh - alpha[0];
  const double qnM = huM * nx + hvM * ny;
  const double qnP = huP * nx + hvP * ny;
  const double qvM = -huM * ny + hvM * nx;
  const double qvP = -huP * ny + hvP * nx;
  double unM, vnM, unP, vnP;
  evaluateVelocity(hcrit, hM, qnM, qvM, &unM, &vnM);
  evaluateVelocity(hcrit, hP, qnP, qvP, &unP, &vnP);
  alpha[0] = 0.5 * (hP - hM - roe->h / roe->c * (unP - unM));
  alpha[1] = roe->h * (vnP - vnM);
  alpha[2] = 0.5 * (hP - hM + roe->h / roe->c * (unP - unM));
#ifdef DEBUG
  mexPrintf("Wave strength\n");
  mexPrintf("local velocity [%f, %f]\n", uM, vM);
  mexPrintf("neigh velocity [%f, %f]\n", uP, vP);
  mexPrintf("alpha = [%f, %f, %f]\n", alpha[0], alpha[1], alpha[2]);
#endif
  return;
}

void evaluateRoeSolver(const double hmin,   ///< water depth threshold
                       const double gra,    ///< gravity acceleration
                       const double hM,     ///< local water depth
                       const double huM,    ///< local flux
                       const double hvM,    ///< local flux
                       const double hP,     ///< neighbour water depth
                       const double huP,    ///< neighbour flux
                       const double hvP,    ///< neighbour flux
                       const double nx,     ///< outward normal vector
                       const double ny,     ///< outward normal vector
                       const RoeState *roe, ///< roe averaged states
                       double *Fh,          ///< roe flux on h
                       double *Fhu,         ///< roe flux on hu
                       double *Fhv          ///< roe flux on hv
) {
  double EM[3], GM[3];
  evaluateFluxTerm2d(hmin, gra, hM, huM, hvM, EM, GM);
  Fh[0] = EM[0] * nx + GM[0] * ny;
  Fhu[0] = EM[1] * nx + GM[1] * ny;
  Fhv[0] = EM[2] * nx + GM[2] * ny;
  evaluateFluxTerm2d(hmin, gra, hP, huP, hvP, EM, GM);
  Fh[0] += EM[0] * nx + GM[0] * ny;
  Fhu[0] += EM[1] * nx + GM[1] * ny;
  Fhv[0] += EM[2] * nx + GM[2] * ny;

  double alpha[3];
  evaluateRoeWaveStrength(hmin, hM, huM, hvM, hP, huP, hvP, nx, ny, roe, alpha);
  const double unroe = roe->u * nx + roe->v * ny;
  const double lambda1 = fabs(unroe - roe->c);
  const double lambda2 = fabs(unroe);
  const double lambda3 = fabs(unroe + roe->c);

#ifdef DEBUG
  mexPrintf("eigenvalue lambda = [%f, %f, %f]\n", lambda1, lambda2, lambda3);
#endif
  Fh[0] -= lambda1 * alpha[0];
  Fhu[0] -= lambda1 * alpha[0] * (roe->u - roe->c * nx);
  Fhv[0] -= lambda1 * alpha[0] * (roe->v - roe->c * ny);

  Fhu[0] += lambda2 * alpha[1] * ny;
  Fhv[0] -= lambda2 * alpha[1] * nx;

  Fh[0] -= lambda3 * alpha[2];
  Fhu[0] -= lambda3 * alpha[2] * (roe->u + roe->c * nx);
  Fhv[0] -= lambda3 * alpha[2] * (roe->v + roe->c * ny);

  Fh[0] *= 0.5;
  Fhu[0] *= 0.5;
  Fhv[0] *= 0.5;

#ifdef DEBUG
  mexPrintf("Roe flux = [%f, %f, %f]\n", Fh[0], Fhu[0], Fhv[0]);
#endif
  return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  FluxSolver solver = ConvertInputMexVariable2d(nlhs, nrhs, plhs, prhs);

  const size_t NdimOut = 3;
  const size_t K = solver.K;
  const size_t TNfp = solver.TNfp;
  const mwSize dimOut[3] = {TNfp, K, 3};
  plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double *Fh = mxGetPr(plhs[0]);
  double *Fqx = Fh + TNfp * K;
  double *Fqy = Fh + 2 * TNfp * K;

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
      const double hvM = solver.hvM[sk];
      const double hvP = solver.hvP[sk];
#ifdef DEBUG
      mexPrintf("k = %d, sk = %d\n", k, n);
      mexPrintf("h = [%f, %f]\nhu = [%f, %f]\nhv = [%f, %f]\n", hM, hP, huM,
                huP, hvM, hvP);
#endif
      if ((hM > solver.hmin) || (hP > solver.hmin)) {
        const double nx = solver.nx[sk];
        const double ny = solver.ny[sk];
        RoeState roe;
        evaluateRoeAverage(solver.gra, solver.hmin, hM, huM, hvM, hP, huP, hvP,
                           &roe);
        evaluateRoeSolver(solver.hmin, solver.gra, hM, huM, hvM, hP, huP, hvP,
                          nx, ny, &roe, Fh + sk, Fqx + sk, Fqy + sk);
      }
    }
  }
  return;
}