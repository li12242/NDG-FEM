/**
 * \file Evaluate Lax–Friedrichs numerical flux function
 * \details
 * For more details, please see Giraldo et al. (2002).
 *
 * [1] Giraldo FX, Hesthaven JS, Warburton T. Nodal High-Order Discontinuous
 * Galerkin Methods for the Spherical Shallow Water Equations. Journal of
 * Computational Physics 2002;181:499–525.
 */

#include "../../SWEAbstractNumFluxSolver2d/private/SWENumFlux2d.h"

void evaluateLocalJacobian(
    const double hmin, ///< wet depth threshold
    const double gra,  ///< gravity acceleration
    const double h,    ///< depth
    const double hu,   ///< flux
    const double hv,   ///< flux
    double *Eh,        ///< flux term
    double *Ehu,       ///< flux term
    double *Ehv,       ///< flux term
    double *lambda     ///< maximum eigenvalues of flux Jacobian
) {
  double u, v;
  if (h > hmin) {
    u = hu / h;
    v = hv / h;
    // evaluate local eigenvalues of the flux Jacobian
    *lambda = sqrt(u * u + v * v) + sqrt(gra * h);
  } else {
    u = 0.0;
    v = 0.0;
    *lambda = 0.0;
  }
  const double huv = h * u * v;
  const double h2 = h * h;
  *Eh = hu;
  *Ehu = h * u * u + 0.5 * gra * h2;
  *Ehv = huv;
  return;
}

void evaluateLFFlux(const double hmin, ///< water threshold
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
  double lambM, lambP, EM[3], EP[3];
  evaluateLocalJacobian(hmin, gra, hM, qnM, qvM, EM + 0, EM + 1, EM + 2,
                        &lambM);
  evaluateLocalJacobian(hmin, gra, hP, qnP, qvP, EP + 0, EP + 1, EP + 2,
                        &lambP);
  double lambda = max(lambM, lambP);
  *Fhn = 0.5 * (EM[0] + EP[0] - lambda * (hP - hM));
  *Fqxn = 0.5 * (EM[1] + EP[1] - lambda * (qnP - qnM));
  *Fqyn = 0.5 * (EM[2] + EP[2] - lambda * (qvP - qvM));
  return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  FluxSolver solver = ConvertInputMexVariable2d(nlhs, nrhs, plhs, prhs);

  const size_t NdimOut = 3;
  const mwSize dimOut[3] = {solver.TNfp, solver.K, 3};
  plhs[0] = mxCreateNumericArray(NdimOut, dimOut, mxDOUBLE_CLASS, mxREAL);

  double *Fh = mxGetPr(plhs[0]);
  double *Fqx = Fh + solver.TNfp * solver.K;
  double *Fqy = Fh + 2 * solver.TNfp * solver.K;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
  for (int k = 0; k < solver.K; k++) {
    for (int n = 0; n < solver.TNfp; n++) {
      const size_t sk = k * solver.TNfp + n;

      const double nx = solver.nx[sk];
      const double ny = solver.ny[sk];

      double qnM, qnP, qvM, qvP;
      RotateFluxToNormal2d(solver.huM[sk], solver.hvM[sk], nx, ny, &qnM, &qvM);
      RotateFluxToNormal2d(solver.huP[sk], solver.hvP[sk], nx, ny, &qnP, &qvP);
      double Fqns, Fqvs;
      evaluateLFFlux(solver.hmin, solver.gra, solver.hM[sk], qnM, qvM,
                     solver.hP[sk], qnP, qvP, Fh + sk, &Fqns, &Fqvs);

      RotateNormalFluxToCoordinate2d(Fqns, Fqvs, nx, ny, Fqx + sk, Fqy + sk);
    }
  }
  return;
}