#include "mxSWE2d.h"
#include "mex.h"
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void evaluateExactAdjacentNodeValue(NdgEdgeType type,
                                    const double nx, const double ny,
                                    double *fm, double *fp, double *fe){
    switch (type){
        case NdgEdgeInner:
            fe[0] = fp[0];
            fe[1] = fp[1];
            fe[2] = fp[2];
            break;
        case NdgEdgeSlipWall:
            evaluateSlipWallAdjacentNodeValue(nx, ny, fm, fe);
            break;
        case NdgEdgeNonSlipWall:
            evaluateNonSlipWallAdjacentNodeValue(nx, ny, fm, fe);
            break;
        case NdgEdgeZeroGrad:
            fe[0] = fm[0]; fe[1] = fm[1]; fe[2] = fm[2];
            break;
        case NdgEdgeClamped:
            fe[0] = 2*fe[0] - fm[0];
            fe[1] = 2*fe[1] - fm[1];
            fe[2] = 2*fe[2] - fm[2];
            break;
        case NdgEdgeClampedDepth:
            fe[0] = 2*fe[0] - fm[0];
            fe[1] = fm[1];
            fe[2] = fm[2];
            break;
        case NdgEdgeClampedVel:
            fe[0] = fm[0];
            fe[1] = 2*fe[1] - fm[1];
            fe[2] = 2*fe[2] - fm[2];
            break;
        case NdgEdgeFlather:
            evaluateFlatherAdjacentNodeValue(nx, ny, fm, fe);
            break;
        default:
            break;
    }
    return;
}

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

/**
 * @brief Calculation of the HLL numerical flux
 * @param [in] hmin - threadhold of the water depth
 * @param [in] gra - gravity accelerated
 * @param [in] hM, qnM, qvM - water depth on local element node
 * @param [in] hP, qnP, qvP - water depth on adjacent element node
 * @param [out] Fhn, Fqxn, Fqyn - HLL numerical flux;
 */
void evaluateHLLFormula(double hmin,
                        double gra,
                        double hM,
                        double qnM,
                        double qvM,
                        double hP,
                        double qnP,
                        double qvP,
                        double *Fhn,
                        double *Fqxn,
                        double *Fqyn) {
    
    double Em[3], Ep[3];
    double Gm[3], Gp[3];
    evaluateFluxTerm(hmin, gra, hM, qnM, qvM, Em, Gm);
    evaluateFluxTerm(hmin, gra, hP, qnP, qvP, Ep, Gp);
//    reduce_nodal_flux(hmin, gra, hM, qnM, qvM, &EhM, &EqnM, &EqvM, &GhM, &GqnM,
//                      &GqvM);
//    reduce_nodal_flux(hmin, gra, hP, qnP, qvP, &EhP, &EqnP, &EqvP, &GhP, &GqnP,
//                      &GqvP);
    /* calculation of wave speed */
    double sM, sP, us, cs, unM, uvM, unP, uvP;
//    evaluateFlowRateByDeptheThreshold(hmin, hM, qnM, qvM, &unM, &uvM);
//    evaluateFlowRateByDeptheThreshold(hmin, hP, qnP, qvP, &unP, &uvP);
    if ((hM > hmin) & (hP > hmin)) {
        unM = qnM / hM;
        unP = qnP / hP;
        us = (double)(0.5 * (unM + unP) + sqrt(gra * hM) - sqrt(gra * hP));
        cs = (double)(0.5 * (sqrt(gra * hM) + sqrt(gra * hP)) + 0.25 * (unM - unP));
        
        sM = (double)min(unM - sqrt(gra * hM), us - cs);
        sP = (double)max(unP + sqrt(gra * hP), us + cs);
    } else if ((hM > hmin) & (hP <= hmin)) {
        unM = qnM / hM;
        sM = (double)(unM - sqrt(gra * hM));
        sP = (double)(unM + 2 * sqrt(gra * hM));
    } else if ((hM <= hmin) & (hP > hmin)) {
        unP = qnP / hP;
        sM = (double)(unP - 2 * sqrt(gra * hP));
        sP = (double)(unP + sqrt(gra * hP));
    } else { /* both dry element */
        sM = 0;
        sP = 0;
    }
    /* HLL flux function */
    if ((sM >= 0) & (sP > 0)) {
        *Fhn = Em[0];
        *Fqxn = Em[1];
        *Fqyn = Em[2];
    } else if ((sM < 0) & (sP > 0)) {
        *Fhn = (sP * Em[0] - sM * Ep[0] + sM * sP * (hP - hM)) / (sP - sM);
        *Fqxn = (sP * Em[1] - sM * Ep[1] + sM * sP * (qnP - qnM)) / (sP - sM);
        *Fqyn = (sP * Em[2] - sM * Ep[2] + sM * sP * (qvP - qvM)) / (sP - sM);
    } else if ((sM < 0) & (sP <= 0)) {
        *Fhn = Ep[0];
        *Fqxn = Ep[1];
        *Fqyn = Ep[2];
    } else if ( (fabs(sM) < EPS) & (fabs(sP) < EPS) ) {
        *Fhn = Em[0];
        *Fqxn = Em[1];
        *Fqyn = Em[2];
    } else {
        mexErrMsgIdAndTxt("Matlab:mxEvaluateFlux2d:ErrWaveSpeed",
                          "The wave speed computation occurs an error.");
    }
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* check input & output */
    if (nrhs != 9){
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
    double *eidP = mxGetPr(prhs[3]);
    signed char *eidtype = (signed char *)mxGetData(prhs[4]);
    double *nx = mxGetPr(prhs[5]);
    double *ny = mxGetPr(prhs[6]);
    double *fext = mxGetPr(prhs[7]);
    double *fphys = mxGetPr(prhs[8]);
    
    const mwSize *dims = mxGetDimensions(prhs[8]);
    
    const size_t Np = dims[0];
    const size_t K = dims[1];
    const size_t TNfp = mxGetM(prhs[2]);
    
    const size_t ndimOut = 3;
    const mwSize dimOut[3] = {TNfp, K, Nvar};
    plhs[0] = mxCreateNumericArray(ndimOut, dimOut, mxDOUBLE_CLASS, mxREAL);
    
    double *dFh = mxGetPr(plhs[0]);
    double *dFqx = dFh + TNfp*K;
    double *dFqy = dFh + 2*TNfp*K;
    
    double *h = fphys;
    double *hu = fphys + K*Np;
    double *hv = fphys + 2*K*Np;
    double *z = fphys + 3*K*Np;
    
    double *he = fext;
    double *hue = fext + K*Np;
    double *hve = fext + 2*K*Np;
    double *ze = fext + 3*K*Np;
    
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
            const int ip = (int) eidP[sk] - 1;
            
            double fm[4] = {h[im], hu[im], hv[im], z[im]};
            double fp[4] = {h[ip], hu[ip], hv[ip], z[ip]};
            double fe[4] = {he[im], hue[im], hve[im], ze[im]};
            
            const double nx_ = nx[sk];
            const double ny_ = ny[sk];
            
            evaluateExactAdjacentNodeValue(type, nx_, ny_, fm, fp, fe);
            
            const double qnM =  fm[1] * nx_ + fm[2] * ny_;
            const double qvM = -fm[1] * ny_ + fm[2] * nx_;
            const double qnP =  fe[1] * nx_ + fe[2] * ny_;
            const double qvP = -fe[1] * ny_ + fe[2] * nx_;
            
            double Fhns, Fqns, Fqvs;
            evaluateHLLFormula(hcrit, gra, fm[0], qnM, qvM, fe[0], qnP, qvP, &Fhns, &Fqns, &Fqvs);
            
            dFh[sk]  = -Fhns;
            dFqx[sk] = -(Fqns * nx_ - Fqvs * ny_);
            dFqy[sk] = -(Fqns * ny_ + Fqvs * nx_);
            
            double E[3], G[3];
            evaluateFluxTerm(hcrit, gra, fm[0], fm[1], fm[2], E, G);
            
            dFh[sk] += nx_ * E[0] + ny_ * G[0];
            dFqx[sk] += nx_ * E[1] + ny_ * G[1];
            dFqy[sk] += nx_ * E[2] + ny_ * G[2];
        }
    }
}