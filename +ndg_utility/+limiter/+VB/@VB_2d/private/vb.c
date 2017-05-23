#include "vb.h"
#define DEBUG 0

/**
 * @brief get local variable field from local gradient
 *
 * @param [in] Np number of points
 * @param [in] hmean mean value
 * @param [in] xc,yc centre coordinate
 * @param [in] x,y coordinate
 * @param [in] phpx,phpy element gradient
 * @return h variable value on each nodes
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


void VertLimit(int K, int Np, int Nfaces, int Nfp, double *fmax, double *fmin, 
	double *fc, double *xc, double *yc, 
	double *f_Q, double *x, double *y,
	double *Fmask, double *EToV, double *flim, Wei_Fun_t WeiGrad){

	double *fv = (double*) calloc(Nfaces, sizeof(double));
	double *xv = (double*) calloc(Nfaces, sizeof(double));
	double *yv = (double*) calloc(Nfaces, sizeof(double));
	
	int k,f,n;
	double dfdx, dfdy;
	for (k=0;k<K;k++){
		double xm = xc[k];
		double ym = yc[k];
		double fm = fc[k];
		int tflag = 0;
		for(f=0;f<Nfaces;f++){
			int n = k*Np + (int)Fmask[f*Nfp]-1; // node index
			int v = (int)EToV[k*Nfaces + f]-1; // vertex index
			xv[f] = x[n];
			yv[f] = y[n];
			fv[f] = f_Q[n];

			#if DEBUG
			mexPrintf("k=%d, f=%d, xv=%f, yv=%f, fv=%f, fmax=%f, fmin=%f, ",
				k, f, xv[f], yv[f], fv[f], fmax[v], fmin[v]);
			#endif
			if(fv[f]>fmax[v]){ 
				fv[f]=fmax[v]; tflag = 1; 
			}else if(fv[f]<fmin[v]){
				fv[f]=fmin[v]; tflag = 1; 
			}
			#if DEBUG
			mexPrintf("fv_=%f\n", fv[f]);
			#endif
		}
		if (tflag){
			GetWeiGrad(Nfaces, xv, yv, fv, xm, ym, fm, &dfdx, &dfdy, WeiGrad);
			Grad2Node(Np, fm, xm, ym, x+k*Np, y+k*Np, dfdx, dfdy, flim+k*Np);
		}else{
			for(n=0;n<Np;n++){
				int sk = k*Np +n;
				flim[sk] = f_Q[sk];
			}
		}
		
	}

	free(fv);
	free(xv);
	free(yv);
	return;
}

void GetWeiGrad(int Nsub, double *xv, double *yv, double *hv,
    double xc, double yc, double hc, 
    double *dhdx, double *dhdy, Wei_Fun_t WeiGrad){

    double *gra_x, *gra_y, *gra_det;
    gra_x = (double *) calloc(Nsub, sizeof(double));
    gra_y = (double *) calloc(Nsub, sizeof(double));
    gra_det = (double *) calloc(Nsub, sizeof(double));
    
    #if DEBUG
    mexPrintf("xc=%f, yc=%f, hc=%f\n", xc, yc, hc);
    #endif
    double a[4],x[2],f[2];
    // double frac = Nsub*eps;
    int i;
    for(i=0;i<Nsub;i++){
        /* vertex index */
        int l1 = i;
        int l2 = (i+1)%Nsub;
        /* coefficient matrix and rhs */
        a[0] = xv[l1] - xc; a[1] = yv[l1] - yc;
        a[2] = xv[l2] - xc; a[3] = yv[l2] - yc;
        f[0] = hv[l1] - hc; f[1] = hv[l2] - hc;
        
        #if DEBUG
        mexPrintf("Nfaces=%d, x1=%f, y1=%f, h1=%f, x2=%f, y2=%f, h2=%f\n",
            i, xv[l1], yv[l1], hv[l1], xv[l2], yv[l2], hv[l2]);
        #endif
        /* get local gradient x=(dhdx, dhdy) of ith subdomain */
        MatrixSolver2(a, f, x);
        gra_x[i] = x[0]; gra_y[i] = x[1];
        gra_det[i] = x[0]*x[0] + x[1]*x[1];
        
        #if DEBUG
        mexPrintf("Nfaces=%d, dhdx=%f, dhdy=%f, det=%f\n",
            i, gra_x[i], gra_y[i], gra_det[i]);
        #endif
    }

    WeiGrad(Nsub, gra_x, gra_y, gra_det, dhdx, dhdy);

    free(gra_x);
    free(gra_y);
    free(gra_det);
    return;
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
