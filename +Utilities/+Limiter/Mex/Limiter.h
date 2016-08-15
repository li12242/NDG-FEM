#include "mex.h"
#include <math.h>

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

#define real double

void minmod(int n, real *a, real *m);
real sign  (real a);