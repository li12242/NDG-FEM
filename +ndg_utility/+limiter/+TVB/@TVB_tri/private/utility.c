#include "TVB.h"

/**
* @brief solve for equations with 2 unknows
* @details solve the equation of \f[A \cdot x = f \f], while the coefficient
* matrix A is
* \f[ A = \begin{bmatrix} a[0], & a[1] \cr
* a[2], & a[3] \end{bamtrix} \f].
* The equations is solved by multiply the inverse matrix
* \f[A^{-1} = \frac{1}{\left\| A \right\|}\begin{bmatrix} a[3], & -a[1] \cr
* -a[2], & a[0] \end{bamtrix}\f] to the rhs vector f, giving by
* \f[ x=A^{-1} \cdot f \f], while \f[ \left\| A \right\| = a[0]a[3] - a[1]a[2] \f$]
* is the norm of matrix.
*
* @param [double] a coefficient matrix
* @param [double] f rhs vector
* @return
* name     | type     | description of value
* -------- |----------|----------------------
* x   | double      | result for equations
*
*/
void MatrixSolver2(double *a, double *f, double *x){

   double det = a[0]*a[3] - a[1]*a[2];
   x[0] = ( f[0]*a[3] - f[1]*a[1])/det;
   x[1] = (-f[0]*a[2] + f[1]*a[0])/det;
   return;
}

/**
 * @brief get local variable field from local gradient
 *
 * @param [int] Np number of points
 * @param [double] hmean mean value
 * @param [double] xc centre coordinate
 * @param [double] yc centre coordinate
 * @param [double*] x coordinate
 * @param [double*] y coordinate
 * @param [double] phpx element gradient
 * @param [double] phpy element gradient
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * h   | double | variable value on each nodes
 *
 */
void Grad2Node(int Np, double hmean,
    double xc, double yc,
 	double *x, double *y,
    double phpx, double phpy, double *h)
{
 	int i;
 	for(i=0;i<Np;i++){
            double dx = x[i] - xc;
            double dy = y[i] - yc;
            h[i] = hmean + dx*phpx + dy*phpy;
    }
}

/**
 * @brief Minmod function
 *
 * @details
 * The minmod function
 * \f[m(a_1, a_2, \cdots, a_m) = \left\{ \begin{array}{ll}
 * s \mathrm{min}_{1\le i\le m} \left| a_i \right|, & \left|s \right| = 1 \cr
 * 0, & otherwise,
 * \end{array} \right. \quad s = \frac{1}{m}\sum_{i=1}^m{sign(a_i)} \f]
 *
 * @param[in] n number of vector
 * @param[in] a vector
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * m   | double* | the limiter result
 *
 */
void minmod(int n, double *a, double *m){
	int i;
	double f, s, t;
	/* use the first element to initialize */
	t = fabs(a[0]);
	f = sign(a[0]);
	s = f;

	for(i=1;i<n;i++){
		f = sign(a[i]);

		if (s*f<0.0){
			*m = 0.0;
			return;
		}
		t = min( t, fabs(a[i]) );
	}

	*m = t*f;
	return;
}
