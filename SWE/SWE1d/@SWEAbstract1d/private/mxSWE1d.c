#include "mxSWE1d.h"

void evaluateSurfFluxTerm(double hmin,
    double gra,
    double h,
    double hu,
    double* E)
{
    double u;
    evaluateFlowRateByDeptheThreshold(hmin, h, hu, &u);
    const double h2 = h * h;
    E[0] = hu;
    E[1] = h * u * u + 0.5 * gra * h2;
    return;
}

void evaluateFlowRateByDeptheThreshold(double hcrit,
    double h,
    double hu,
    double* u)
{
    if (h > hcrit) {
        *u = hu / h;
    } else {
        *u = 0;
        //    double h4 = h * h * h * h;
        //    *u = sqrt(2.0) * h * hu / sqrt(h4 + max(h4, TAU));
    }
    return;
}

void evaluateFlowRateByCellState(const NdgRegionType type,
    const double h,
    const double hu,
    double* u)
{
    switch (type) {
    case NdgRegionWet:
        *u = hu / h;
        break;
    default: {
        *u = 0.0;
        //      double h4 = h * h * h * h;
        //      *u = sqrt(2.0) * h * hu / sqrt(h4 + max(h4, TAU));
        break;
    }
    }
}

void evaluateSlipWallAdjacentNodeValue(const double nx,
    double* fm,
    double* fp)
{
    const double hM = fm[0];
    const double qxM = fm[1];
    double qnM = qxM * nx; // outward normal flux
    // adjacent value
    fp[0] = hM;
    fp[1] = (-qnM) * nx;
    return;
}

void evaluateNonSlipWallAdjacentNodeValue(const double nx,
    double* fm,
    double* fp)
{
    const double hM = fm[0];
    const double qxM = fm[1];
    // adjacent value
    fp[0] = hM;
    fp[1] = -qxM;
    return;
}

void evaluateFlatherAdjacentNodeValue(double nx, double* fm, double* fe)
{
    const double gra = 9.81;
    double hM = fm[0];
    double h_ext = fe[0];
    double qx_ext = fe[1];
    double qn_ext = qx_ext * nx; // outward normal flux

    double qn = qn_ext - sqrt(gra * hM) * (h_ext - hM);
    // double qn = - sqrt(solver.gra*hM)*(2*etaE - etaM);
    fe[0] = fm[0];
    fe[1] = qn * nx;
    // f_P[1] = qn*nx;
    // f_P[2] = qn*ny;
    return;
}