#ifndef VECTOROPERATOR_H
#define VECTOROPERATOR_H

#include "mex.h"

typedef struct{
	double x, y; // coordinate
}POINT;

double cross(POINT A, POINT B);
double dot(POINT A, POINT B);
void minus(POINT A, POINT B, POINT *C);

#endif