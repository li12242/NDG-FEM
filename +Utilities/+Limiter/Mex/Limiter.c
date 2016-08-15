#include "Limiter.h"
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
 * m   | real* |
 *
 */
void minmod(int n, real *a, real *m){
	int i;
	real f, s, t;
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

/**
 * @brief sign function
 *
 * @param[in] a real number
 *
 * @return    
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * s   | real |
 */
real sign(real a){
	real m;
	if (a>0.0){
		return 1.0;
	}else if (a<0.0){
		return -1.0;
	}else{
		return 0.0;
	}
}