#include "mxSWE2d.h"
#include "mex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/** */
void evaluateNodalFlux(double hcrit, double gra,
                double h, double qx, double qy, double z,
                double *Eh, double *Eqx, double *Eqy,
                double *Gh, double *Gqx, double *Gqy)
{
    if (h > hcrit) {
        double h2 = 0.5 * gra * (h * h - z * z);
        double huv = qx * qy / h;
        
        *Eh = qx;
        *Gh = qy;
        *Eqx = (qx * qx / h + h2);
        //*Eqx = qx*qx/h;
        *Gqx = huv;
        *Eqy = huv;
        *Gqy = (qy * qy / h + h2);
        //*Gqy = qy*qy/h;
    } else { // for dry nodes
        *Eh = 0;
        *Eqx = 0;
        *Eqy = 0;
        *Gh = 0;
        *Gqx = 0;
        *Gqy = 0;
    }
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* check input & output */
    if (nrhs != 4){
        mexErrMsgIdAndTxt("Matlab:mxEvaluateFlux2d:InvalidNumberInput",
                          "4 inputs required.");
    }
    
    if (nlhs != 2){
        mexErrMsgIdAndTxt("Matlab:mxEvaluateFlux2d:InvalidNumberOutput",
                          "2 output required.");
    }
    
    double hcrit = mxGetScalar(prhs[0]);
    double gra = mxGetScalar(prhs[1]);
    signed char *regType = (signed char *)mxGetData(prhs[2]);
    double *fphys = mxGetPr(prhs[3]);
    
    const mwSize *dims = mxGetDimensions(prhs[3]);
    
    const size_t Np = dims[0];
    const size_t K = dims[1];
    
    const size_t ndimOut = 3;
    const mwSize dimOut[3] = {Np, K, Nvar};
    plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
    
    double *E = mxGetPr(plhs[0]);
    double *G = mxGetPr(plhs[1]);
    
    double *Eh = E;
    double *Ehu = E + K*Np;
    double *Ehv = E + 2*K*Np;
    double *Gh = G;
    double *Ghu = G + K*Np;
    double *Ghv = G + 2*K*Np;

    
    double *h = fphys;
    double *hu = fphys + K*Np;
    double *hv = fphys + 2*K*Np;
    double *z = fphys + 3*K*Np;

#ifdef _OPENMP
#pragma omp parallel for num_threads(DG_THREADS)
#endif
    for (int k=0; k<K; k++) {
        NdgRegionType type = (NdgRegionType) regType[k];
        if (type == NdgRegionDry)
            continue;
        
        for (int n=0; n<Np; n++) {
            size_t sk = k*Np + n;
            evaluateNodalFlux(hcrit, gra, h[sk], hu[sk], hv[sk], z[k],
                       Eh+sk, Ehu+sk, Ehv+sk, Gh+sk, Ghu+sk, Ghv+sk);

//            double h_ = h[sk];
//            double hu_ = hu[sk];
//            double hv_ = hv[sk];
//            double z_ = z[sk];
            
            //            double h2 = h_*h_ - z_*z_;
//            double u_, v_;
//            evaluateFlowRateByCellState(type, h_, hu_, hv_, &u_, &v_);
//            
//            double huv = h_*u_*v_;
//            Eh[sk] = hu_;
//            Gh[sk] = hv_;
//            Ehu[sk] = h_*u_*u_ + 0.5*gra*h2;
//            Ghu[sk] = huv;
//            Ehv[sk] = huv;
//            Ghv[sk] = h_*v_*v_ + 0.5*gra*h2;
        }
    }
    
    return;
}