#include "mxSWE2d.h"
#include <math.h>

/**
 Evaluate the flow rate depending on the depth threshold
 */
void evaluateFlowRateByDeptheThreshold(const double hcrit,
                                       const double h, const double hu, const double hv,
                                       double *u, double *v){
    
    if(h>hcrit){
        *u = hu/h;
        *v = hv/h;
    }else{
        *u = 0.0;
        *v = 0.0;
    }

    return;
}

/**
 Evaluate the flow rate depending on the cell states
 */
void evaluateFlowRateByCellState(const NdgRegionType type,
                                 const double h, const double hu, const double hv,
                                 double *u, double *v){
    switch (type) {
        case NdgRegionWet:
            *u = hu/h; *v = hv/h;
            break;
        default:
            *u = 0; *v = 0;
            break;
    }
}


void evaluateSlipWallAdjacentNodeValue(const double nx, const double ny, double *fm, double *fp){
    const double hM  = fm[0];
    const double qxM = fm[1];
    const double qyM = fm[2];
    double qnM =  qxM*nx + qyM*ny; // outward normal flux
    double qvM = -qxM*ny + qyM*nx; // outward tangential flux
    // adjacent value
    fp[0] = hM;
    fp[1] = (-qnM)*nx - qvM*ny;
    fp[2] = (-qnM)*ny + qvM*nx;
    return;
}

void evaluateNonSlipWallAdjacentNodeValue(const double nx, const double ny, double *fm, double *fp){
    const double hM  = fm[0];
    const double qxM = fm[1];
    const double qyM = fm[2];
    // adjacent value
    fp[0] = hM;
    fp[1] = -qxM;
    fp[2] = -qyM;
    return;
}

void evaluateFlatherAdjacentNodeValue(double nx, double ny, double *fm, double *fe){
    const double gra = 9.81;
    double hM = fm[0];
    double h_ext = fe[0];
    double qx_ext = fe[1];
    double qy_ext = fe[2];
    double qn_ext =  qx_ext*nx + qy_ext*ny; // outward normal flux
    double qv_ext = -qx_ext*ny + qy_ext*nx; // tangential flux
    
    double qn = qn_ext - sqrt(gra*hM)*(h_ext - hM);
    double qv = qv_ext;
    // double qn = - sqrt(solver.gra*hM)*(2*etaE - etaM);
    fe[0] = fm[0];
    fe[1] = qn*nx - qv*ny;
    fe[2] = qn*ny + qv*nx;
    // f_P[1] = qn*nx;
    // f_P[2] = qn*ny;
    return;
}