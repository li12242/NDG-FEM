/**
 * @file Limiter.c
 * @brief
 * The model for NDG-FEM limiter functions
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */
#include "Limiter.h"

/**
 * @brief Calculate the mean gradient through n sub-domian
 * @param [int] Nsub number of subdomain
 * @param [real*] xv coordinate of vertex
 * @param [real*] xv coordinate of vertex
 * @param [real*] hv
 * @param [real] xc
 * @param [real] yc
 * @param [real] hcs
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * dhdx | real* | averaged gradient for x coordinate
 * dhdx | real* | averaged gradient for y coordinate
 */
void meanGradient(int Nsub, real *xv, real *yv, real *hv,
    real xc, real yc, real hc, real *dhdx, real *dhdy){

    real eps = 1.0e-12;
    real *gra_x, *gra_y, *gra_det;
    gra_x = (real *) calloc(Nsub, sizeof(real));
    gra_y = (real *) calloc(Nsub, sizeof(real));
    gra_det = (real *) calloc(Nsub, sizeof(real));

    real a[4],x[2],f[2];
    real frac = Nsub*eps;
    int i,j;
    for(i=0;i<Nsub;i++){
        /* vertex index */
        int l1 = i;
        int l2 = (i+1)%Nsub;
        /* coefficient matrix and rhs */
        a[0] = xv[l1] - xc; a[1] = yv[l1] - yc;
        a[2] = xv[l2] - xc; a[3] = yv[l2] - yc;
        f[0] = hv[l1] - hc; f[1] = hv[l2] - hc;
        // mexPrintf("Nfaces=%d, x1=%f, y1=%f, h1=%f, x2=%f, y2=%f, h2=%f\n",
        //     i, xv[l1], yv[l1], hv[l1], xv[l2], yv[l2], hv[l2]);
        /* get local gradient x=(dhdx, dhdy) of ith subdomain */
        matrixSolver2(a, f, x);
        gra_x[i] = x[0]; gra_y[i] = x[1];
        gra_det[i] = x[0]*x[0] + x[1]*x[1];
        // mexPrintf("Nfaces=%d, dhdx=%f, dhdy=%f, det=%f\n",
        //     i, gra_x[i], gra_y[i], gra_det[i]);

        // frac += pow(gra_det[i], (Nsub-1.0));
    }
    // mexPrintf("frac=%f\n", frac);

    *dhdx = 0.0; *dhdy = 0.0;
    for(i=0;i<Nsub;i++){
        real w = 1.0;
        for(j=0;j<Nsub;j++){
            if(j==i) {continue;}
            w *= gra_det[j];
        }
        frac += w;
        w += eps;
        // mexPrintf("Nfaces=%d, w=%f\n", i, w);
        *dhdx += w*gra_x[i];
        *dhdy += w*gra_y[i];
    }
    *dhdx /= frac;
    *dhdy /= frac;

    free(gra_x);
    free(gra_y);
    free(gra_det);
    return;
}

/**
 * @brief get local variable field distribution from local gradient
 *
 * @param [int] Np number of points
 * @param [real] hmean mean value
 * @param [real] xc centre coordinate
 * @param [real] yc centre coordinate
 * @param [real*] x coordinate
 * @param [real*] y coordinate
 * @param [real] phpx element gradient
 * @param [real] phpy element gradient
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * h   | real | variable value on each nodes
 *
 */
void getLocalVar(int Np, real hmean, real xc, real yc,
 	real *x, real *y, real phpx, real phpy, real *h){
 	int i;
 	for(i=0;i<Np;i++){
            real dx = x[i] - xc;
            real dy = y[i] - yc;

            h[i] = hmean + dx*phpx + dy*phpy;
    }
}

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
 * @param [real] a coefficient matrix
 * @param [real] f rhs vector
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * x   | real      | result for equations
 *
 */
void matrixSolver2(real *a, real *f, real *x){

	real det = a[0]*a[3] - a[1]*a[2];
	x[0] = ( f[0]*a[3] - f[1]*a[1])/det;
	x[1] = (-f[0]*a[2] + f[1]*a[0])/det;

	return;
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
 * @param[in] Nfaces number of faces
 * @param[in] Nfp number of nodes on each face
 * @param[in] Nfields number of variables for calculation
 * @param[in] h variables in element
 * @param[in] ws integral coefficient
 * @param[in] sJ Jacobian of spicific ith faces
 * @param[in] Fmask local node index of spicific ith faces
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * face_mean   | real[Nfields] | average value on face
 * face_len    | real | length of face
 *
 * Usages:
 * 	 for(f=0;f<Nfaces;f++){
 *	 	face_mean(Nfaces, Nfp, h+k*Np, ws, sJ+k*Nfaces*Nfp+f*Nfp, Fmask+f,
 *			&face_mean, &face_len);
 *   }
 *
 */
void faceMean(int Nfaces, int Nfp, int Nfields, real **h,
	real *ws, real *sJ, real *Fmask,
	real *face_mean, real *face_len){
	/* initialization */
	*face_len  = 0.0;
    int fld;
    for(fld=0;fld<Nfields;fld++){
        face_mean[fld] = 0.0;
    }

    int fnp,sp;
	for(fnp=0;fnp<Nfp;fnp++){
	  	sp = (int) (*(Fmask + fnp*Nfaces) - 1); // local index of face
        real w=ws[fnp], j=sJ[fnp];
        // mexPrintf("fnp=%d, w*j=%f, ", fnp, w*j);
        for(fld=0;fld<Nfields;fld++){
            // mexPrintf("var[%d]=%f, fmean[%d]=%f, ", fld, *(h[fld]+sp), fld, face_mean[fld]);
            face_mean[fld] += w*j*( *(h[fld]+sp));
            // mexPrintf("var[%d]=%f, fmean[%d]=%f, ", fld, *(h[fld]+sp), fld, face_mean[fld]);
        }
        // mexPrintf("\n");
		*face_len += w*j;
	}
    for(fld=0;fld<Nfields;fld++){
        face_mean[fld] /= *face_len;
    }
	return;
}

/**
 * @brief calculation of cell averages
 *
 * @param[int] Np number of points in each cell
 * @param[int] Nfields number of variables for calculation
 * @param[real**] h variables in element
 * @param[real*] w integral coefficient
 * @param[real*] J Jacobian of spicific ith faces
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * cellMean | real[Nfields] | average value of each variable field
 * area | real* | area
 *
 * Usages:
 * 	 for(k=0;k<K;k++){
 *     fld[0] = h+k*Np; fld[1] = x+k*Np; fld[2] = y+k*Np;
 *     cellMean(Np, 3, fld, w, J+k*Np, cmean, area+k);
 *     hmean[k] = cmean[0];
 *     xmean[k] = cmean[1];
 *     ymean[k] = cmean[2];
 *   }
 *
 */
void cellMean(int Np, int Nfields, real **h,
        real *w, real *J, real *cellMean, real *area){

    int i,fld;
    //initialize
    *area  = 0.0;
    for(fld=0;fld<Nfields;fld++){
        cellMean[fld] = 0.0;
    }
    //loop over nodes
    for(i=0;i<Np;i++){
        real wi = w[i];
        real Ji = J[i];
        for(fld=0;fld<Nfields;fld++){
            cellMean[fld] += wi*Ji*( *(h[fld]+i));
        }
        *area  += wi*Ji;
    }
    for(fld=0;fld<Nfields;fld++){
        cellMean[fld] /= *area;
    }
    return;
}
