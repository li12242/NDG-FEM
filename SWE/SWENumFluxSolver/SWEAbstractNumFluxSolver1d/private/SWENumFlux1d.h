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
  double *hM;
  double *hP;
  double *huM;
  double *huP;
} FluxSolver;

#define NRHS 5

/** Put input variable into FluxSolver  */
FluxSolver ConvertInputMexVariable1d(const int Nlhs,       ///< number of LHS
                                     const int Nrhs,       ///< number of RHS
                                     const mxArray *prhs[] ///< RHS pointer
) {
  /* check input & output */
  if (Nrhs != NRHS) {
    mexPrintf("Matlab:%s:InvalidNumberInput,\n", __FILE__);
    mexPrintf("%d inputs required.\n", NRHS);
  }

  FluxSolver solver;
  solver.hmin = mxGetScalar(prhs[0]);
  solver.gra = mxGetScalar(prhs[1]);
  solver.nx = mxGetPr(prhs[2]);
  solver.hM = mxGetPr(prhs[3]);
  solver.hP = mxGetPr(prhs[4]);

  const mwSize *dims = mxGetDimensions(prhs[4]);
  const size_t TNfp = dims[0];
  const size_t K = dims[1];
  const size_t NFP = K * TNfp;
  solver.TNfp = TNfp;
  solver.K = K;
  solver.huM = solver.hM + NFP;
  solver.huP = solver.hP + NFP;

  return solver;
}

/** Evaluate flux term in surface integration */
void evaluateFluxTerm1d(const double hmin, ///< water threshold
                        const double gra,  ///< gravity acceleration
                        const double h,    ///< water depth
                        const double hu,   ///< flux variable
                        double *E          ///< surface integral flux term
) {
  double u;
  if (h > hmin) {
    u = hu / h;
  } else {
    u = 0.0;
  }
  const double h2 = h * h;
  E[0] = hu;
  E[1] = h * u * u + 0.5 * gra * h2;
  return;
}

/** Rotate flux to outward normal direction */
void RotateFluxToNormal1d(const double hu, ///< flux at x component
                          const double nx, ///< outward normal vector
                          double *qn       ///< normal flux
) {
  *qn = +hu * nx;
  return;
}

/** Rotate normal flux to coordinate direction */
void RotateNormalFluxToCoordinate1d(const double Fqn, ///< normal flux
                                    const double nx,  ///< outward normal vector
                                    double *Fqx       ///< flux at x component
) {
  *Fqx = Fqn * nx;
  return;
}

#endif // SWENUMFLUX2D_H