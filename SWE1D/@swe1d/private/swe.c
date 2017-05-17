#include "swe.h"

/*
 *
 */
void nodal_flux(double hcrit, double gra,
    double h, double qx, double *Fh, double *Fq)
{
    if(h>hcrit){
        *Fh = qx;
        *Fq = (qx*qx/h + 0.5*gra*h*h);
    }else{
        *Fh = 0;
        *Fq = 0;
    }
    return;
}
