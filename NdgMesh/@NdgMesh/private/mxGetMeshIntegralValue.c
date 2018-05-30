/*
 * @file
 * Calculate the integral averaged value in each cells.
 *
 * @details
 * Calculate the cell intergal averaged values by the formula
 * \f$ \int_{\Omega_k} u dA = \sum_{n=1}^{N_q} w_n J_n u(\xi_n) \f$
 * where \f$ \xi_n \f$ is the Gauss quadrature points.
 *
 * @param[in] fvar The field value on each interpolation points in all elements
 * @param[in] wq The quadrature weights in each cell
 * @param[in] J The Jacobian determination on interpolation points in all
 * elements
 * @param[in] Vq The transform matrix to convert nodal values to quadrature
 * values
 * @return fint The integral values in each elements
 *
 * @code
 *  [ fint ] = mxGetMeshIntegralValue( fvar, wq, J, Vq )
 * @endcode
 *
 */
#include "mex.h"
#include <math.h>
#include "blas.h"

#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* check input & output */
    if (nrhs != 4)
        mexErrMsgIdAndTxt("Matlab:mxGetIntegralValue:InvalidNumberInput",
                          "4 inputs required.");
    if (nlhs != 1)
        mexErrMsgIdAndTxt("Matlab:mxGetIntegralValue:InvalidNumberOutput",
                          "1 output required.");
    
    /* get inputs */
    double *fvar = mxGetPr(prhs[0]);
    double *wq = mxGetPr(prhs[1]);
    double *J = mxGetPr(prhs[2]);
    double *Vq = mxGetPr(prhs[3]);
    
    /* get dimensions */
    size_t Np = mxGetM(prhs[0]); // number of interpolation points
    size_t Nq = mxGetM(prhs[3]); // number of quadrature points
    size_t K = mxGetN(prhs[0]);  // number of elements
    
    /* allocate output array */
    plhs[0] = mxCreateDoubleMatrix((mwSize)1, (mwSize)K, mxREAL);
    double *fint = mxGetPr(plhs[0]);
    
    const double one = 1;
    const double zero = 0;
    ptrdiff_t one_ptrdiff = 1;
    ptrdiff_t Np_ptrdiff = Np;
    ptrdiff_t Nq_ptrdiff = Nq;
    char *tran = "N";
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    
    for (int k = 0; k < K; k++) {
        double Jq[Nq], fq[Nq];
        
        // map the node values fvar to quadrature nodes by
        // \f$ fq = Vq * fvar \f$
        dgemm(tran, tran, &Nq_ptrdiff, &one_ptrdiff, &Np_ptrdiff, &one, Vq,
              &Nq_ptrdiff, fvar + k * Np, &Np_ptrdiff, &zero, fq, &Nq_ptrdiff);
        dgemm(tran, tran, &Nq_ptrdiff, &one_ptrdiff, &Np_ptrdiff, &one, Vq,
              &Nq_ptrdiff, J + k * Np, &Np_ptrdiff, &zero, Jq, &Nq_ptrdiff);
        
        for (int n = 0; n < Nq; n++) {
            fint[k] += wq[n] * Jq[n] * fq[n];
        }
    }
}
