#include "mex.h"
#include "blas.h"
#include "conv2d.h"

#define DEBUG 0

void inner_fv_term(int Np, int K, double *h, double *u, double *v,
                   int Nedge, double *v1, double *v2,
                   double *nx, double *ny, double *ds,
                   signed char *EToR, double *rhs)
{
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k = 0; k < K; k++)
    {
        if ((cell_type)EToR[k] != REFINE)
            continue;

        int n, ind = k * Nedge;
        for (n = 0; n < Nedge; n++)
        {
            int n1 = (int)v1[n] + k*Np - 1; // change to C type
            int n2 = (int)v2[n] + k*Np - 1;
            double nx_ = nx[ind];
            double ny_ = ny[ind];
            double delta_s = ds[ind];

            double num_flux;
            upwind_flux(h[n1], h[n2], u[n1], v[n1], nx_, ny_, &num_flux);

            rhs[n1] -= num_flux * delta_s;
            rhs[n2] += num_flux * delta_s;
#if DEBUG
            mexPrintf("k=%d, n=%d, n1=%d, n2=%d, nx=%e, ny=%e, delta=%e, flux=%e, rhs1=%e, rhs2=%e\n",
                      k, n, n1, n2, nx_, ny_, delta_s, num_flux, rhs[n1], rhs[n2]);
#endif
            ind++;
        }
    }

#if DEBUG
    mexPrintf("inner_rhs = \n");
    int n;
    for (n = 0; n < Np; n++)
    {
        mexPrintf("\t");
        for (k = 0; k < K; k++)
        {
            mexPrintf("%e\t", rhs[k * Np + n]);
        }
        mexPrintf("\n");
    }
#endif
    return;
}

void surface_fv_term(size_t Np, size_t Nfp, size_t K,
                     double *h, double *h_ext, double *u, double *v,
                     signed char *EToR, signed char *eidtype,
                     double *eidM, double *eidP,
                     double *nx, double *ny, double *Js, double *LIFT,
                     double *rhs)
{

    int k;
    double *flux = calloc(Nfp * K, sizeof(double));
    // #ifdef _OPENMP
    // #pragma omp parallel for num_threads(DG_THREADS)
    // #endif
    for (k = 0; k < K; k++)
    {
        if ((cell_type)EToR[k] != REFINE)
            continue;

        int j, ind = k * Nfp;
        for (j = 0; j < Nfp; j++)
        {
            int iM = (int)eidM[ind] - 1; // change index to C type
            int iP = (int)eidP[ind] - 1;
            double f_M = h[iM]; // local and adjacent node values
            double hP = h[iP];
            double uM = u[iM], vM = v[iM];
            // outward normal vector of local element
            double nx_ = nx[ind];
            double ny_ = ny[ind];
            double f_ext; // external values on local nodes
            f_ext = h_ext[iM];
            bc_type type = (bc_type)eidtype[ind];
            // get adjacent values hP, qxP, qyP, considering
            // various boudnary conditions
            double f_P;
            int info = bound_cond(f_M, hP, f_ext, nx_, ny_, type, &f_P);
            // if(info) mexErrMsgTxt("Unknown boundary conditions.");
            double numflux;
            upwind_flux(f_M, f_P, uM, vM, nx_, ny_, &numflux);

            flux[ind] = -numflux;
            ind++;
        }
    }
#if DEBUG
    mexPrintf("flux = \n");
    int n;
    for (n = 0; n < Nfp; n++)
    {
        mexPrintf("\t");
        for (k = 0; k < K; k++)
        {
            mexPrintf("%e\t", flux[k * Nfp + n]);
        }
        mexPrintf("\n");
    }
#endif
    double *stemp = calloc(Nfp * K, sizeof(double));
    dvecm(Nfp * K, 1, flux, Js, stemp);
    char *chn = "N";
    double one = 1.0;
    mwSignedIndex Np_ = Np, K_ = K, Nfp_ = Nfp;
    dgemm(chn, chn, &Np_, &K_, &Nfp_, &one, LIFT, &Np_, stemp, &Nfp_, &one, rhs, &Np_);
    free(flux);
    free(stemp);
    return;
}

/*
 * @brief calculate the sub-cell discretization by finite volume scheme.
 * Usages: 
 *      [ rhsQ ] = rhs_fv_term(f_Q, f_ext, u, v, Nedge, v1, v2, nx, ny, ds, ...
 *                              eidM, eidP, nxM, nyM, Js, LIFT, vol, ...
 *                              eidtype, EToR)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* check input */
    if (nrhs != 19)
    {
        mexErrMsgIdAndTxt("MATLAB:rhs_fv_term:invalidNumInputs",
                          "20 input required.");
    }
    else if (nlhs > 1)
    {
        mexErrMsgIdAndTxt("MATLAB:rhs_fv_term:maxlhs",
                          "Too many output arguments.");
    }

    double *f_Q = mxGetPr(prhs[0]);
    double *f_extQ = mxGetPr(prhs[1]);
    double *u = mxGetPr(prhs[2]);
    double *v = mxGetPr(prhs[3]);
    double Nedge = mxGetScalar(prhs[4]);
    double *v1 = mxGetPr(prhs[5]);
    double *v2 = mxGetPr(prhs[6]);
    double *nx = mxGetPr(prhs[7]);
    double *ny = mxGetPr(prhs[8]);
    double *ds = mxGetPr(prhs[9]);
    double *eidM = mxGetPr(prhs[10]);
    double *eidP = mxGetPr(prhs[11]);
    double *nxM = mxGetPr(prhs[12]);
    double *nyM = mxGetPr(prhs[13]);
    double *Js = mxGetPr(prhs[14]);
    double *LIFT = mxGetPr(prhs[15]);
    double *vol = mxGetPr(prhs[16]);
    signed char *eidtype = (signed char *)mxGetData(prhs[17]);
    signed char *EToR = (signed char *)mxGetData(prhs[18]);

    /* get dimensions */
    size_t Np = mxGetM(prhs[0]);
    size_t K = mxGetN(prhs[0]);
    size_t Nfp = mxGetM(prhs[10]);

    /* allocate output array */
    plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
    double *rhsQ = mxGetPr(plhs[0]);

    double *vtemp = calloc(Np * K, sizeof(double));
    inner_fv_term(Np, K, f_Q, u, v,
                  (int)Nedge, v1, v2, nx, ny, ds, EToR, vtemp);

    surface_fv_term(Np, Nfp, K, f_Q, f_extQ, u, v,
                    EToR, eidtype, eidM, eidP, nxM, nyM, Js, LIFT, vtemp);

    dvecd(Np * K, 1, vtemp, vol, rhsQ);

    free(vtemp);
    return;
}