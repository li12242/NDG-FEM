#include "pbswe.h"

/*
 *
 */
void nodal_flux(double hcrit, double gra,
    double eta, double qx, double bot, double *Fh, double *Fq)
{
    double h = eta - bot;
    if(h>hcrit){
        *Fh = qx;
        *Fq = ( qx*qx/h + 0.5*gra*(eta*eta - 2*eta*bot) );
    }else{
        *Fh = 0;
        *Fq = 0;
    }
    return;
}
