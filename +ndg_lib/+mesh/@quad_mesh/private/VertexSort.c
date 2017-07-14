#include "VertexSort.h"

POINT stdvert;

// Cross product function
double Multiply(POINT p1, POINT p2, POINT p3) {
	return ((p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x));
}

// Distance of points
double Distance(POINT p1, POINT p2){
	return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
}

// Compare function
int cmp(const void *p1, const void *p2){
	POINT *p3, *p4;
	double m;
	p3 = (POINT *)p1;
	p4 = (POINT *)p2;
	m = Multiply(stdvert, *p3, *p4);
	if (m < 0) return 1;
	else if (m == 0 && (Distance(stdvert, *p3) < Distance(stdvert, *p4)))
		return 1;
	else return -1;
}


void VertexSort(int Nvert, double *EToV, 
        double *x, double *y, double *newEToV){
	int i;
	POINT vertex[Nvert];

	for(i=0;i<Nvert;i++){
		int t = (int) EToV[i] - 1;
		vertex[i].x = x[t];
		vertex[i].y = y[t];
		vertex[i].nv = t;
	}
	stdvert.x  = vertex[0].x;
	stdvert.y  = vertex[0].y;
	stdvert.nv = vertex[0].nv;

	qsort(&vertex[1], Nvert-1, sizeof(POINT), cmp);

	for(i=0;i<Nvert;i++){
		newEToV[i] = (double) vertex[i].nv + 1;
	}
	return;
}