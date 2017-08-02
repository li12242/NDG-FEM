#include "conv2d.h"

int nodal_flux(double c, double u, double v, double *E, double *G){
	*E = u*c; 
	*G = v*c;
	return 0;
}