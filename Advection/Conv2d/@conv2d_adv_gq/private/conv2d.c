#include "conv2d.h"

int nodal_flux(double c, double u, double v, double *E, double *G) {
  *E = u * c;
  *G = v * c;
  return 0;
}

int bound_cond(double varM, double varP, double f_ext, double nx, double ny,
               bc_type type, double *f_P) {
  switch (type) {
  case Inner:
  case SlipWall:
  case NSlipWall:
    *f_P = varP;
    break;
  case ZeroGrad:
    *f_P = varM;
    break;
  case Clamped:
  case ClampedDepth:
  case ClampedVel:
  case Flather:
    *f_P = 2 * f_ext - varM;
    break;
  }

  return 0;
}