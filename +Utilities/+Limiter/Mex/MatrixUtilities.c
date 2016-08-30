#include "Limiter.h"

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
 * @param [int] a coefficient matrix
 * @param [int] f rhs vector
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * x   | real      | result for equations
 * 
 */
void MatrixSolver2(real *a, real *f, real *x){

	real det = a[0]*a[3] - a[1]*a[2];
	x[0] = ( f[0]*a[3] - f[1]*a[1])/det;
	x[1] = (-f[0]*a[2] + f[1]*a[0])/det;

	return;
}