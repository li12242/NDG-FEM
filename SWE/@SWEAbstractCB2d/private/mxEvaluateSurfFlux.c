#include "mxSWE2d.h"
#include "mex.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void evaluateFluxTerm(const double hmin, const double gra,
                      const double h, const double hu, const double hv,
                      double *E, double *G){
    
    double u, v;
    evaluateFlowRateByDeptheThreshold(hmin, h, hu, hv, &u, &v);
    const double huv = h*u*v;
    const double h2 = h*h;
    E[0] = hu;
    G[0] = hv;
    E[1] = h*u*u + 0.5*gra*h2;
    G[1] = huv;
    E[2] = huv;
    G[2] = h*v*v + 0.5*gra*h2;
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* check input & output */
    if (nrhs != 7){
        mexErrMsgIdAndTxt("Matlab:mxEvaluateFlux2d:InvalidNumberInput",
                          "4 inputs required.");
    }
    
    if (nlhs != 1){
        mexErrMsgIdAndTxt("Matlab:mxEvaluateFlux2d:InvalidNumberOutput",
                          "2 output required.");
    }
    
    double hcrit = mxGetScalar(prhs[0]);
    double gra = mxGetScalar(prhs[1]);
    double *eidM = mxGetPr(prhs[2]);
    signed char *eidtype = (signed char *)mxGetData(prhs[3]);
    double *nx = mxGetPr(prhs[4]);
    double *ny = mxGetPr(prhs[5]);
    double *fphys = mxGetPr(prhs[6]);
    
    const mwSize *dims = mxGetDimensions(prhs[6]);
    
    const size_t Np = dims[0];
    const size_t K = dims[1];
    const size_t TNfp = mxGetM(prhs[2]);
    
    const size_t ndimOut = 3;
    const mwSize dimOut[3] = {TNfp, K, Nvar};
    plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
    
    double *Fhs = mxGetPr(plhs[0]);
    double *Fqxs = Fhs + TNfp*K;
    double *Fqys = Fhs + 2*TNfp*K;
    
    double *h = fphys;
    double *hu = fphys + K*Np;
    double *hv = fphys + 2*K*Np;
    double *z = fphys + 3*K*Np;
    
#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k<K; k++) {
        for (int n=0; n<TNfp; n++) {
            const size_t sk = k*TNfp + n;
            NdgEdgeType type = (NdgEdgeType)eidtype[sk];
            if ( type == NdgEdgeGaussEdge ) {
                continue;
            }
            
            const int im = (int) eidM[sk] - 1;
            double fm[4] = {h[im], hu[im], hv[im], z[im]};
            
            const double nx_ = nx[sk];
            const double ny_ = ny[sk];
            
            double E[3], G[3];
            evaluateFluxTerm(hcrit, gra, fm[0], fm[1], fm[2], E, G);
            
            Fhs[sk] = nx_ * E[0] + ny_ * G[0];
            Fqxs[sk] = nx_ * E[1] + ny_ * G[1];
            Fqys[sk] = nx_ * E[2] + ny_ * G[2];
        }
    }
}