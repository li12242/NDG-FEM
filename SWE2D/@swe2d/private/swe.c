#include "swe.h"

/*
 *
 */
void nodal_flux(double hcrit, double gra,
    double h, double qx, double qy, double z,
    double *Eh, double *Eqx, double *Eqy,
    double *Gh, double *Gqx, double *Gqy)
{
    if(h>hcrit){
        double h2 = 0.5*gra*(h*h - z*z);
        double huv = qx*qy/h;

        *Eh = qx;
        *Gh = qy;
        *Eqx = (qx*qx/h + h2);
        *Gqx = huv;
        *Eqy = huv;
        *Gqy = (qy*qy/h + h2);
    }else{ // for dry nodes
        *Eh = 0; *Eqx = 0; *Eqy = 0;
        *Gh = 0; *Gqx = 0; *Gqy = 0;
    }
    return;
}
