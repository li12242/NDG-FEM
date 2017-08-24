#include "nconv2d.h"

int nodal_flux(double c, double *E, double *G)
{
    *E = c * c * 0.5;
    *G = c * c * 0.5;
    return 0;
}

int bound_cond(double varM, double varP, double f_ext,
               double nx, double ny, bc_type type,
               double *f_P)
{

    switch (type)
    {
    case Inner:
    case SlipWall:
    case NSlipWall:
        *f_P = varP;
        break;
    case ZeroGrad:
        *f_P = varM;
        break;
    case Clamped:
        *f_P = f_ext;
        break;
    case ClampedDepth:
    case ClampedVel:
    case Flather:
        *f_P = 2 * f_ext - varM;
        break;
    }
    return 0;
}