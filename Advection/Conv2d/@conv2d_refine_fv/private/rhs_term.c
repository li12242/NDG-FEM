#include "mex.h"
#include "blas.h"
#include "conv2d.h"

#define DEBUG 0

/*
 */
void rhs_parall(size_t Np, size_t K, size_t Nfp,
                double *Dr, double *Ds, double *LIFT,
                double *rx, double *ry, double *sx, double *sy,
                double *J, double *Js,
                double *dflux, double *eflux, double *gflux,
                double *rhs)
{
    int k;
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (k = 0; k < K; k++)
    {
        /* volume term */
        int n, m;
        for (m = 0; m < Np; m++)
        {
            double rx_ = rx[k * Np + m];
            double ry_ = ry[k * Np + m];
            double sx_ = sx[k * Np + m];
            double sy_ = sy[k * Np + m];

            double *rhsQ = rhs + (k * Np + m);
            double *E = eflux + k * Np;
            double *G = gflux + k * Np;
            for (n = 0; n < Np; n++)
            {
                double dr = Dr[n * Np + m];
                double ds = Ds[n * Np + m];
                double dx = rx_ * dr + sx_ * ds;
                double dy = ry_ * dr + sy_ * ds;
                rhsQ[0] -= dx * E[n] + dy * G[n];
            }

            double j = 1 / J[k * Np + m];
            for (n = 0; n < Nfp; n++)
            {
                double js = Js[k * Nfp + n];
                double dfs = dflux[k * Nfp + n];

                rhsQ[0] += LIFT[m + n * Np] * dfs * js * j;
            }
        }
    }
    return;
}
/*
 */
void rhs_term(size_t Np_, size_t K_, size_t Nfp_,
              double *Dr, double *Ds, double *LIFT,
              double *rx, double *ry, double *sx, double *sy,
              double *J, double *Js,
              double *dflux, double *eflux, double *gflux,
              double *rhs)
{
    char *chn = "N";
    double one = 1.0, zero = 0.0;
    double *vtemp = calloc(Np_ * K_, sizeof(double));
    mwSignedIndex Np = Np_, K = K_, Nfp = Nfp_;
    dgemm(chn, chn, &Np, &K, &Np, &one, Dr, &Np, eflux, &Np, &zero, vtemp, &Np);
    dvecm(Np * K, -1, rx, vtemp, rhs);
    dgemm(chn, chn, &Np, &K, &Np, &one, Ds, &Np, eflux, &Np, &zero, vtemp, &Np);
    dvecm(Np * K, -1, sx, vtemp, rhs);
    dgemm(chn, chn, &Np, &K, &Np, &one, Dr, &Np, gflux, &Np, &zero, vtemp, &Np);
    dvecm(Np * K, -1, ry, vtemp, rhs);
    dgemm(chn, chn, &Np, &K, &Np, &one, Ds, &Np, gflux, &Np, &zero, vtemp, &Np);
    dvecm(Np * K, -1, sy, vtemp, rhs);

    double *stemp = calloc(Nfp * K, sizeof(double));
    dvecm(Nfp * K, 1, Js, dflux, stemp);
    dgemm(chn, chn, &Np, &K, &Nfp, &one, LIFT, &Np, stemp, &Nfp, &zero, vtemp, &Np);
    dvecd(Np * K, 1, vtemp, J, rhs);

    free(vtemp);
    free(stemp);
    return;
}

/**
 * @brief calculate the R.H.S of conv2d problem.
 *
 * Usages:
 *  [ rhsQ ] = rhs_term(f_Q, f_ext, u, v, 
 *                      nx, ny, eidM, eidP, eidtype, EToR,    % for surface term
 *                      Dr, Ds, rx, ry, sx, sy, LIFT, J, Js)  % for rhs term
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* check input */
    if (nrhs != 19)
    {
        mexErrMsgIdAndTxt("MATLAB:rhs_term:invalidNumInputs",
                          "20 input required.");
    }
    else if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("MATLAB:rhs_term:maxlhs",
                          "Too many output arguments.");
    }

    double *f_Q = mxGetPr(prhs[0]);
    double *f_extQ = mxGetPr(prhs[1]);
    double *u = mxGetPr(prhs[2]);
    double *v = mxGetPr(prhs[3]);
    double *nx = mxGetPr(prhs[4]);
    double *ny = mxGetPr(prhs[5]);
    double *eidM = mxGetPr(prhs[6]);
    double *eidP = mxGetPr(prhs[7]);
    signed char *eidtype = (signed char *)mxGetData(prhs[8]);
    signed char *EToR = (signed char *)mxGetData(prhs[9]);
    double *Dr = mxGetPr(prhs[10]);
    double *Ds = mxGetPr(prhs[11]);
    double *rx = mxGetPr(prhs[12]);
    double *ry = mxGetPr(prhs[13]);
    double *sx = mxGetPr(prhs[14]);
    double *sy = mxGetPr(prhs[15]);
    double *LIFT = mxGetPr(prhs[16]);
    double *J = mxGetPr(prhs[17]);
    double *Js = mxGetPr(prhs[18]);

    /* get dimensions */
    size_t Np = mxGetM(prhs[0]);
    size_t K = mxGetN(prhs[0]);
    size_t Nfp = mxGetM(prhs[6]);

    /* allocate output array */
    plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
    double *rhsQ = mxGetPr(plhs[0]);

    /* surfce integral term */
    double *dflux = calloc(Nfp * K, sizeof(double));
    surf_term(Nfp, K, f_Q, f_extQ, u, v, nx, ny, eidM, eidP, eidtype, EToR, dflux);

    // #if DEBUG
    //     mexPrintf("dflux = \n");
    //     int n, k;
    //     for (n = 0; n < Nfp; n++)
    //     {
    //         mexPrintf("\t");
    //         for (k = 0; k < K; k++)
    //         {
    //             mexPrintf("%e\t", dflux[k * Np + n]);
    //         }
    //         mexPrintf("\n");
    //     }
    // #endif
    /* volume flux term */
    double *eflux = calloc(Np * K, sizeof(double));
    double *gflux = calloc(Np * K, sizeof(double));
    flux_term(Np, K, f_Q, u, v, EToR, eflux, gflux);

    // #if DEBUG
    //     mexPrintf("e = \n");
    //     for (n = 0; n < Np; n++)
    //     {
    //         mexPrintf("\t");
    //         for (k = 0; k < K; k++)
    //         {
    //             mexPrintf("%f\t", eflux[k * Np + n]);
    //         }
    //         mexPrintf("\n");
    //     }

    //     mexPrintf("g = \n");
    //     for (n = 0; n < Np; n++)
    //     {
    //         mexPrintf("\t");
    //         for (k = 0; k < K; k++)
    //         {
    //             mexPrintf("%f\t", gflux[k * Np + n]);
    //         }
    //         mexPrintf("\n");
    //     }
    // #endif
    /* DG rhs term */
    rhs_term(Np, K, Nfp, Dr, Ds, LIFT, rx, ry, sx, sy, J, Js, dflux, eflux, gflux, rhsQ);
    /* fv surface term */
    free(dflux);
    free(eflux);
    free(gflux);
    return;
}