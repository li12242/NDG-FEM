//
//  mxVertLimit
//
//  Created by li12242 on 17/10/31.
//  Copyright (c) 2017年 li12242. All rights reserved.
//

#include "mex.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define EPSILON 1.0e-12

void evaluateWenoLocalGrad(size_t Nsub,
    double* subGfx,
    double* subGfy,
    double* subGraDet,
    double* gfx,
    double* gfy)
{
    double frac = 0.0;
    double r = 2.0; // a positive number
    *gfx = 0.0;
    *gfy = 0.0;
    for (int i = 0; i < Nsub; i++) {
        double w = pow(sqrt(subGraDet[i]) + EPSILON, -r);
        frac += w;
        *gfx += w * subGfx[i];
        *gfy += w * subGfy[i];
        // if(k==29 | k==149)
        //     mexPrintf("k=%d, w[%d]=%f\n", k, i, w);
    }
    *gfx /= frac;
    *gfy /= frac;
    // if(k==29 | k==149)
    //     mexPrintf("k=%d, gx=%f, gy=%f\n", k, *gfx, *gfy);
    // return;
}

/* the weights of van Albada limiter */
void evaluateVALocalGrad(
	int Nsub, 
	double *gra_x, 
	double *gra_y, 
	double *gra_det,
    double *dhdx, 
    double *dhdy)
{
    double frac=Nsub*EPSILON;;
    int i,j;

    for(*dhdx=0.0,*dhdy=0.0,i=0;i<Nsub;i++){
        double w = 1.0;
        for(j=0;j<Nsub;j++){
            if(i==j) continue;
            w = w*gra_det[j];
        }
        w += EPSILON;
        frac += w;
        *dhdx += w*gra_x[i];
        *dhdy += w*gra_y[i];
    }
    *dhdx /= frac;
    *dhdy /= frac;
}

/* weights of Hermite WENO limiter */
void evaluateJKLocalGrad(int Nsub, double *gra_x, double *gra_y, double *gra_det,
    double *dhdx, double *dhdy){
    double frac=Nsub*EPSILON;
    int i,j;
    for(i=0;i<Nsub;i++){ frac += (pow(gra_det[i], (Nsub-1.0)) + EPSILON ); }

    for(*dhdx=0.0,*dhdy=0.0,i=0;i<Nsub;i++){
        double w = 1.0;
        for(j=0;j<Nsub;j++){
            if(i==j) continue;
            w = w*gra_det[j];
        }
        w += EPSILON;
        *dhdx += w*gra_x[i];
        *dhdy += w*gra_y[i];
    }
    *dhdx /= frac;
    *dhdy /= frac;
}

/**
 * @brief 
 * Solve for equations with 2 unknows.
 *
 * @details 
 * Solve the equation of \f[A \cdot x = f \f], 
 * while the coefficient matrix A is
 * \f[ A = \begin{bmatrix} a[0], & a[1] \cr a[2], & a[3] \end{bamtrix} \f].
 *
 * The equations is solved by multiply the inverse matrix
 * \f[A^{-1} = \frac{1}{\left\| A \right\|}\begin{bmatrix} a[3], & -a[1] \cr
 * -a[2], & a[0] \end{bamtrix}\f] 
 * to the rhs vector f, giving by
 * \f[ x=A^{-1} \cdot f \f], while \f[ \left\| A \right\| = a[0]a[3] - a[1]a[2] \f$]
 * is the norm of matrix.
 *
 * @param [in] a The coefficient matrix
 * @param [in] f The RHS vector
 * @param [out] x Solutions
 */
void MatrixSolver2(double* a, double* f, double* x)
{

    double det = a[0] * a[3] - a[1] * a[2];
    x[0] = ( f[0] * a[3] - f[1] * a[1]) / det;
    x[1] = (-f[0] * a[2] + f[1] * a[0]) / det;
    return;
}

void evaluateVertexWeightedGradient(size_t Nsub,
    double* cellvx,
    double* cellvy,
    double* cellfv,
    double xc,
    double yc,
    double fc,
    double* gfx,
    double* gfy)
{
    double subGfx[Nsub];
    double subGfy[Nsub];
    double subGraDet[Nsub];
    double a[4], x[2], f[2];
    // double frac = Nsub*eps;
    for (int n = 0; n < Nsub; n++) {
        /* vertex index */
        int l1 = n;
        int l2 = (n + 1) % Nsub;
        /* coefficient matrix and rhs */
        a[0] = cellvx[l1] - xc;
        a[1] = cellvy[l1] - yc;
        a[2] = cellvx[l2] - xc;
        a[3] = cellvy[l2] - yc;
        f[0] = cellfv[l1] - fc;
        f[1] = cellfv[l2] - fc;

        /* get local gradient x=(dhdx, dhdy) of ith subdomain */
        MatrixSolver2(a, f, x);
        subGfx[n] = x[0];
        subGfy[n] = x[1];
        subGraDet[n] = x[0] * x[0] + x[1] * x[1];
    }
    evaluateWenoLocalGrad(Nsub, subGfx, subGfy, subGraDet, gfx, gfy);
        // if (k==29 | k==149){
        //    for( int n = 0; n < Nsub; n++){
        //         mexPrintf("k=%d, subGfx[%d]=%f, subGfy[%d]=%f\n", k, n, subGfx[n], n, subGfy[n]);
        //     }
        // }
    return;
}

/**
 * @brief Get interpolation node values from the gradient and cell averages.
 *
 * @param [in] Np Number of interpolations
 * @param [in] fmean cell integral averaged value
 * @param [in] xc,yc centre coordinate
 * @param [in] x,y coordinate
 * @param [in] gfx,gfy element gradient
 * @param [out] fvar variable value on each nodes
 *
 */
void projectGradToNodeValue(size_t Np,
    double fmean,
    double xc,
    double yc,
    double* x,
    double* y,
    double gfx,
    double gfy,
    double* fvar)
{

    for (int i = 0; i < Np; i++) {
        double dx = x[i] - xc;
        double dy = y[i] - yc;
        fvar[i] = fmean + dx * gfx + dy * gfy;
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

    /* check input & output */
    if (nrhs != 13) {
        mexErrMsgIdAndTxt(
                "Matlab:mxVertLimit:InvalidNumberInput", 
                "8 inputs required.");
    }

    /* get inputs */
    double* fvar = mxGetPr(prhs[0]);
    double* x = mxGetPr(prhs[1]);
    double* y = mxGetPr(prhs[2]);
    double* xc = mxGetPr(prhs[3]);
    double* yc = mxGetPr(prhs[4]);
    double* vx = mxGetPr(prhs[5]);
    double* vy = mxGetPr(prhs[6]);
    double* fvert = mxGetPr(prhs[7]);
    double* fvmin = mxGetPr(prhs[8]);
    double* fvmax = mxGetPr(prhs[9]);
    double* cvar = mxGetPr(prhs[10]);
    double* EToV = mxGetPr(prhs[11]);
    double* Fmask = mxGetPr(prhs[12]);

    /* get dimensions */
    size_t Np = mxGetM(prhs[0]); // number of interpolation points
    size_t Nv = mxGetM(prhs[11]); // number of vertex in each cell
    size_t K = mxGetN(prhs[0]); // number of elements
    size_t Nfp = mxGetM(prhs[12]);

    plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
    double* flimit = mxGetPr(plhs[0]);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k = 0; k < K; k++) {
        double xm = xc[k];
        double ym = yc[k];
        double fm = cvar[k];
        bool troubleCellFlag = 0;
//         bool troubleCellFlag = 1;

        double cellvf[Nv];
        double cellvx[Nv];
        double cellvy[Nv];

        for (int n = 0; n < Nv; n++) {
            size_t nodeId = k * Np + (int)Fmask[n * Nfp] - 1;
            size_t vertId = (int)EToV[k * Nv + n] - 1;
            cellvx[n] = vx[vertId];
            cellvy[n] = vy[vertId];

//             cellvf[n] = fvert[vertId];
            cellvf[n] = fvar[nodeId];

            if (cellvf[n] > fvmax[vertId]) {
                troubleCellFlag = 1;
                cellvf[n] = fvert[vertId];
            } else if (cellvf[n] < fvmin[vertId]) {
                troubleCellFlag = 1;
                cellvf[n] = fvert[vertId];
            }
        }
        if (troubleCellFlag) {
            double gfx, gfy;
            evaluateVertexWeightedGradient(
                Nv, cellvx, cellvy, cellvf, xm, ym, fm, &gfx, &gfy);
            projectGradToNodeValue(
                Np, fm, xm, ym,
                x + k * Np, y + k * Np,
                gfx, gfy, flimit + k * Np);
            // if( k==29 | k==149 )
            //     mexPrintf("k=%d, fm=%f, gfx=%f, gfy=%f\n", k, fm, gfx, gfy);

        } else {
            for (int n = 0; n < Np; n++) {
                flimit[k * Np + n] = fvar[k * Np + n];
            }
        }
    }
    return;
}
