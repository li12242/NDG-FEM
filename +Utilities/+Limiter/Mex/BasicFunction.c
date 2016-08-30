#include "Limiter.h"

void GetLocalVar(int Np, real hmean, real xc, real yc, 
 	real *x, real *y, 
 	real phpx, real phpy, real *h){
 	int i;
 	for(i=0;i<Np;i++){
            real dx = x[i] - xc;
            real dy = y[i] - yc;

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
 * m   | real* | the limiter result
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

/**
 * @brief calculation of averages on faces
 *
 * @param[in] Nfaces No. of faces
 * @param[in] Nfp No. of nodes on each face
 * @param[in] h variables in element
 * @param[in] ws integral coefficient
 * @param[in] sJ Jacobian 
 * @param[in] Fmask 
 *
 * @return    
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * face_mean   | real | average value on face
 * face_len    | real | length of face
 *
 * Usages:
 * 	 for(f=0;f<Nfaces;f++){
 *	 	face_mean(Nfaces, Nfp, h+k*Np, ws, sJ+k*Nfaces*Nfp+f*Nfp, Fmask+f,
 *			&face_mean, &face_len);
 *   }
 *
 */
void FaceMean(int Nfaces, int Nfp, real *h, 
	real *ws, real *sJ, real *Fmask,
	real *face_mean, real *face_len){
	/* initialization */
	*face_len  = 0.0;
    *face_mean = 0.0;

    int fnp,sp;
	for(fnp=0;fnp<Nfp;fnp++){
	  	sp = (int) (*(Fmask + fnp*Nfaces) - 1); // local index of face
		*face_mean    += ws[fnp]*sJ[fnp]*( *(h+sp));
		*face_len  	  += ws[fnp]*sJ[fnp];
	}
	*face_mean /= *face_len;
	return;
}