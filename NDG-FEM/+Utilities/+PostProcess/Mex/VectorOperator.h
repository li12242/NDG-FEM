#ifndef VECTOROPERATOR_H
#define VECTOROPERATOR_H

#include "mex.h"
#include <math.h>

typedef struct{
	double x, y; // coordinate
}POINT;

double cross(POINT A, POINT B);
double dot(POINT A, POINT B);
void minus(POINT A, POINT B, POINT *C);

#endif