#ifndef VERTEXSORT_H
#define VERTEXSORT_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define real double

typedef struct{
	double x, y; // coordinate
	int nv; // No. of point
}POINT;

void VertexSort(int Ne, int Nvert, double *EToV, double *x, double *y, double *newEToV);

#endif