#include "VectorOperator.h"

double cross(POINT A, POINT B){
	double p;
	p = A.x*B.y - B.x*A.y;
	return p;
}

double dot(POINT A, POINT B){
	double p;
	p = A.x*B.x + A.y*B.y;
	return p;
}

void minus(POINT A, POINT B, POINT *C){
	C->x = A.x - B.x;
	C->y = A.y - B.y;
	return;
}