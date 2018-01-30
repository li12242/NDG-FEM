#include "mex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

#define DEBUG 0

void cal_slope_limiter(double fmean, double gx, double gy, double xc, double yc,
                       double *fmask, int Nfp, int Nv,
                       double *x, double *y, int Np,
                       double *fvmax, double *fvmin,
                       double *alpha);

void cal_gg_gradient(double *f_Q, double *x, double *y, int Np,              // node values
                     double *fmask, int Nfp, int Nv, double *ws, double *Js, // edge info
                     double area, double *gx, double *gy);

/*
 * Limit the node values.
 * 
 * Usages:
 *  [ f_Q ] = BJ_limit_2d(f_Q, x, y, ... % node values
 *      c_mean, xc, yc, area, ... % cell values
 *      vmin, vmax, fmask, EToV, 
 *      Js, ws); % vertex info
 *
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* check input & output */
    if (nrhs != 13)
        mexErrMsgTxt("Wrong number of input arguments.");
    if (nlhs != 1)
        mexErrMsgTxt("Wrong number of output arguments");

    /* get inputs */
    double *f_Q = mxGetPr(prhs[0]);
    double *x = mxGetPr(prhs[1]);
    double *y = mxGetPr(prhs[2]);
    double *f_mean = mxGetPr(prhs[3]);
    double *xc = mxGetPr(prhs[4]);
    double *yc = mxGetPr(prhs[5]);
    double *area = mxGetPr(prhs[6]);
    double *v_min = mxGetPr(prhs[7]);
    double *v_max = mxGetPr(prhs[8]);
    double *fmask = mxGetPr(prhs[9]);
    double *EToV = mxGetPr(prhs[10]);
    double *Js = mxGetPr(prhs[11]);
    double *ws = mxGetPr(prhs[12]);

    size_t Np = mxGetM(prhs[0]);  // # of points in each element
    size_t K = mxGetN(prhs[0]);   // # of elements
    size_t Nfp = mxGetM(prhs[9]); // # of points on each edge
    size_t Nv = mxGetN(prhs[9]);  // # of vertex in each element

#if DEBUG
// mexPrintf("Np=%d, K=%d, Nfp=%d, Nv=%d\n", Np, K, Nfp, Nv);
// for (int n = 0; n < Nfp * Nv; n++)
// {
//     mexPrintf("n=%d, Js=%f, ws=%f\n", n, Js[n], ws[n]);
// }
#endif
    /* allocation of output */
    plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
    double *f_limt = mxGetPr(plhs[0]);

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k = 0; k < K; k++)
    {
        int flag = 0;
        double fvmax[Nv], fvmin[Nv]; // the vertex bounds
        for (int v = 0; v < Nv; v++)
        {
            int vid = (int)EToV[k * Nv + v] - 1; // convert to C type
            fvmax[v] = v_max[vid];
            fvmin[v] = v_min[vid];

            vid = (int)fmask[v * Nfp] - 1;
            double fv = f_Q[k * Np + vid];

            if ((fv > fvmax[v]) | (fv < fvmin[v]))
                flag = 1;

#if DEBUG
            mexPrintf("k=%d, v=%d, fv=%f, fvmax=%f, fvmin=%f, flag=%d\n",
                      k, v, fv, fvmax[v], fvmin[v], flag);
#endif
        }

        if (flag == 0)
        {
            for (int n = 0; n < Np; n++)
            {
                int ind = n + k * Np;
                f_limt[ind] = f_Q[ind];
            }
            continue; // non-trobled cell
        }

        double gx, gy; // gradient
        cal_gg_gradient(f_Q + k * Np, x + k * Np, y + k * Np, Np,
                        fmask, Nfp, Nv, ws, Js + k * Nfp * Nv,
                        area[k], &gx, &gy);
#if DEBUG
        mexPrintf("k=%d, dfdx=%f, dfdy=%f\n", k, gx, gy);
#endif
        double alpha;
        cal_slope_limiter(f_mean[k], gx, gy, xc[k], yc[k],
                          fmask, Nfp, Nv,
                          x + k * Np, y + k * Np, Np,
                          fvmax, fvmin, &alpha);
#if DEBUG
        mexPrintf("k=%d, alpha=%f\n", k, alpha);
#endif
        for (int n = 0; n < Np; n++)
        {
            int ind = n + k * Np;
            double xp = x[ind];
            double yp = y[ind];
            f_limt[ind] = f_mean[k] + alpha * (gx * (xp - xc[k]) + gy * (yp - yc[k]));
#if DEBUG
            mexPrintf("f[%d]=%f, ", n, f_limt[ind]);
#endif
        }
#if DEBUG
        mexPrintf("\n");
#endif
    }

    return;
}

/**
 * Calculate the gradient by Green-Gauss formula.
 */
void cal_gg_gradient(double *f_Q, double *x, double *y, int Np,              // node values
                     double *fmask, int Nfp, int Nv, double *ws, double *Js, // edge info
                     double area, double *gx, double *gy)
{
    *gx = 0; // initialize the gradient
    *gy = 0;
    int sk = 0;
    for (int f = 0; f < Nv; f++)
    {
        int v1 = (int)fmask[f * Nfp] - 1;
        int v2 = (int)fmask[((f + 1) % Nv) * Nfp] - 1;

        double dx = x[v2] - x[v1];
        double dy = y[v2] - y[v1];
#if DEBUG
// mexPrintf("f=%d, v1=%d, v2=%d, dx=%f, dy=%f\n", f, v1, v2, dx, dy);
#endif
        for (int n = 0; n < Nfp; n++)
        {
            int node_id = (int)fmask[f * Nfp + n] - 1;
            double j = Js[sk];
            double w = ws[sk++];

            *gx += j * w * dy * f_Q[node_id];
            *gx -= j * w * dx * f_Q[node_id];

#if DEBUG
// mexPrintf("f=%d, n=%d, nid=%d, j=%f, w=%f, gx=%f, gy=%f, area=%f\n",
//           f, n, node_id, j, w, *gx, *gy, area);
#endif
        }
    }

    *gx /= area;
    *gy /= area;
    return;
}

/** 
 * Calculate the slope limiter from the BJ formula.
 * Reference: Kuzmin (2010), Eq. (16).
 */
void cal_slope_limiter(double fmean, double gx, double gy, double xc, double yc,
                       double *fmask, int Nfp, int Nv,
                       double *x, double *y, int Np,
                       double *fvmax, double *fvmin,
                       double *alpha)
{
    *alpha = 1.0; // initialization
    for (int n = 0; n < Nv; n++)
    {
        int ind = (int)fmask[n * Nfp] - 1; // the index of the vertex
        double xv = x[ind];
        double yv = y[ind];
        double fv = fmean + gx * (xv - xc) + gy * (yv - yc);

        if (fv > fvmax[n])
        {
            double temp = min(1, (fvmax[n] - fmean) / (fv - fmean));
            *alpha = min(temp, *alpha);
        }
        else if (fv < fvmin[n])
        {
            double temp = min(1, (fvmin[n] - fmean) / (fv - fmean));
            *alpha = min(temp, *alpha);
        }
#if DEBUG
        mexPrintf("n=%d, fv=%f, fvmax=%f, fvmin=%f, fmean=%f, alpha=%f\n",
                  n, fv, fvmax[n], fvmin[n], fmean, *alpha);
#endif
    }

    return;
}


