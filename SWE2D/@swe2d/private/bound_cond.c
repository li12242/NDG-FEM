#include "swe.h"

void slip_wall(double nx, double ny, double *varM, double *varP);
void non_slip_wall(double nx, double ny, double *varM, double *varP);
void flather(double nx, double ny, double *f_M, double *f_ext, double *f_P);
/*
 * Get adjacent node values with given boundary conditions.
 *
 */
int bound_cond(double *varM, double *varP, double *f_ext,
    double nx, double ny, bc_type type, double *f_P){

    switch (type){
        case Inner:
            f_P[0] = varP[0]; f_P[1] = varP[1]; f_P[2] = varP[2];
            break;
        case SlipWall: slip_wall(nx, ny, varM, f_P); break;
        case NSlipWall: non_slip_wall(nx, ny, varM, f_P); break;
        case ZeroGrad:
            f_P[0] = varM[0]; f_P[1] = varM[1]; f_P[2] = varM[2];
            break;
        case Clamped:
            f_P[0] = 2*f_ext[0] - varM[0]; 
            f_P[1] = 2*f_ext[1] - varM[1]; 
            f_P[2] = 2*f_ext[2] - varM[2];
            break;
        case ClampedDepth:
            f_P[0] = f_ext[0]; f_P[1] = varM[1]; f_P[2] = varM[2];
            break;
        case ClampedVel:
            f_P[0] = varM[0]; f_P[1] = f_ext[1]; f_P[2] = f_ext[2];
            break;
        case Flather: flather(nx, ny, varM, f_ext, f_P); break;
        default:
            return 1;
    }
    return 0;
}

void slip_wall(double nx, double ny, double *varM, double *varP)
{
    const double hM  = varM[0];
    const double qxM = varM[1];
    const double qyM = varM[2];
    double qnM =  qxM*nx + qyM*ny; // outward normal flux
    double qvM = -qxM*ny + qyM*nx; // outward tangential flux
    // adjacent value
    varP[0] = hM;
    varP[1] = (-qnM)*nx - qvM*ny;
    varP[2] = (-qnM)*ny + qvM*nx;
    return;
}

void non_slip_wall(double nx, double ny, double *varM, double *varP){
    const double hM  = varM[0];
    const double qxM = varM[1];
    const double qyM = varM[2];
    // adjacent value
    varP[0] = hM; varP[1] = -qxM; varP[2] = -qyM;
    return;
}

/*
 * Flather boundary condition from Carter and Merrifield (2007).
 *
 */
void flather(double nx, double ny, double *f_M, double *f_ext, double *f_P){
    const double gra = 9.81;
    double hM = f_M[0];
    double h_ext = f_ext[0];
    double qx_ext = f_ext[1];
    double qy_ext = f_ext[2];
    double qn_ext =  qx_ext*nx + qy_ext*ny; // outward normal flux
    double qv_ext = -qx_ext*ny + qy_ext*nx; // tangential flux

    double qn = qn_ext - sqrt(gra*hM)*(h_ext - hM);
    double qv = qv_ext;
    // double qn = - sqrt(solver.gra*hM)*(2*etaE - etaM);
    f_P[0] = f_M[0]; f_P[1] = qn*nx - qv*ny; f_P[2] = qn*ny + qv*nx;
    // f_P[1] = qn*nx;
    // f_P[2] = qn*ny;
    return;
}
