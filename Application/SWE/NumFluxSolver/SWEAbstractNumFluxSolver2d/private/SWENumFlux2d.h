#ifndef SWENUMFLUX2D_H
#define SWENUMFLUX2D_H

#include "mex.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define TOLERR 1e-6
#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

typedef struct {
  int TNfp;
  int K;

  double hmin;
  double gra;
  double *nx;
  double *ny;
  double *hM;
  double *hP;
  double *huM;
  double *huP;
  double *hvM;
  double *hvP;
} FluxSolver;

#define NRHS 6
#define NLHS 1

/** Put input variable into FluxSolver  */
FluxSolver ConvertInputMexVariable2d(const int Nlhs,        ///< number of LHS
                                     const int Nrhs,        ///< number of RHS
                                     const mxArray *plhs[], ///< LHS pointer
                                     const mxArray *prhs[]  ///< RHS pointer
) {
  /* check input & output */
  if (Nrhs != NRHS) {
    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NRHS);
  }

  if (Nlhs != NLHS) {
    mexPrintf("Matlab:%s:InvalidNumberOutput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NLHS);
  }

  FluxSolver solver;
  solver.hmin = mxGetScalar(prhs[0]);
  solver.gra = mxGetScalar(prhs[1]);
  solver.nx = mxGetPr(prhs[2]);
  solver.ny = mxGetPr(prhs[3]);
  solver.hM = mxGetPr(prhs[4]);
  solver.hP = mxGetPr(prhs[5]);

  const mwSize *dims = mxGetDimensions(prhs[4]);
  const size_t TNfp = dims[0];
  const size_t K = dims[1];
  const size_t NFP = K * TNfp;
  solver.TNfp = TNfp;
  solver.K = K;
  solver.huM = solver.hM + NFP;
  solver.hvM = solver.huM + NFP;
  solver.huP = solver.hP + NFP;
  solver.hvP = solver.huP + NFP;

  return solver;
}

/** Evaluate flux term in surface integration */
void evaluateFluxTerm2d(const double hmin, ///< water threshold
                        const double gra,  ///< gravity acceleration
                        const double h,    ///< water depth
                        const double hu,   ///< flux variable
                        const double hv,   ///< flux variable
                        double *E,         ///< surface integral flux term
                        double *G          ///< surface integral flux term
) {
  double u, v;
  if (h > hmin) {
    u = hu / h;
    v = hv / h;
  } else {
    u = 0.0;
    v = 0.0;
  }
  const double huv = h * u * v;
  const double h2 = h * h;
  E[0] = hu;
  G[0] = hv;
  E[1] = h * u * u + 0.5 * gra * h2;
  G[1] = huv;
  E[2] = huv;
  G[2] = h * v * v + 0.5 * gra * h2;
  return;
}

/** Rotate flux to outward normal direction */
void RotateFluxToNormal2d(const double hu, ///< flux at x component
                          const double hv, ///< flux at y component
                          const double nx, ///< outward normal vector
                          const double ny, ///< outward normal vector
                          double *qn,      ///< normal flux
                          double *qv       ///< tangent flux
) {
  *qn = +hu * nx + hv * ny;
  *qv = -hu * ny + hv * nx;
  return;
}

/** Rotate normal flux to coordinate direction */
void RotateNormalFluxToCoordinate2d(const double Fqn, ///< normal flux
                                    const double Fqv, ///< tangent flux
                                    const double nx,  ///< outward normal vector
                                    const double ny,  ///< outward normal vector
                                    double *Fqx,      ///< flux at x component
                                    double *Fqy       ///< flux at y component
) {
  *Fqx = (Fqn * nx - Fqv * ny);
  *Fqy = (Fqn * ny + Fqv * nx);
  return;
}

#endif // SWENUMFLUX2D_H