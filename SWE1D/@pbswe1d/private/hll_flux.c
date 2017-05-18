#include "pbswe.h"

void hll_nodal_flux(double hmin, double gra, double etaM, double etaP,
    double qnM, double qnP, double botM, double botP, double *Fhn, double *Fqxn){

    double FhM, FqnM, FhP, FqnP;
    nodal_flux(hmin, gra, etaM, qnM, botM, &FhM, &FqnM);
    nodal_flux(hmin, gra, etaP, qnP, botP, &FhP, &FqnP);

    double hM = etaM - botM;
    double hP = etaP - botP;
    /* calculation of wave speed */
    double sM, sP;
    if( (hM>hmin) & (hP>hmin) ){
        double unM=qnM/hM;
        double unP=qnP/hP;
        double us = (0.5*(unM + unP)   + sqrt(gra*hM)   - sqrt(gra*hP));
        double cs = (0.5*(sqrt(gra*hM) + sqrt(gra*hP) ) + 0.25*(unM - unP));

        sM = min(unM-sqrt(gra*hM), us-cs);
        sP = max(unP+sqrt(gra*hP), us+cs);
    }else if ( (hM>hmin) & (hP<=hmin) ){
        double unM=qnM/hM;
        sM = (unM - sqrt(gra*hM) );
        sP = (unM + 2*sqrt(gra*hM) );
    }else if ( (hM<=hmin) & (hP>hmin) ){
        double unP=qnP/hP;
        sM = (unP - 2*sqrt(gra*hP) );
        sP = (unP + sqrt(gra*hP) );
    }else{ /* both dry element */
        sM = 0;
        sP = 0;
    }

    /* HLL function */
    if ( (sM>=0) & (sP>0) ){
        *Fhn = FhM; *Fqxn = FqnM;
    }else if((sM<0) & (sP>0)){
        *Fhn  = (sP*FhM  - sM*FhP  + sM*sP*(hP  - hM ))/(sP - sM);
        *Fqxn = (sP*FqnM - sM*FqnP + sM*sP*(qnP - qnM))/(sP - sM);
    }else if( (sM<0)&(sP<=0) ){
        *Fhn = FhP; *Fqxn = FqnP;
    }else if( (sM==0) & (sP==0) ){
        *Fhn = 0; *Fqxn = 0;
    }else{
        mexErrMsgIdAndTxt("SWE_Mex_HLL1d:",
        "The wave speed computation occurs an error.");
    }
    return;
}

/* SWE_Mex_HLL2d
 * Calculation of the HLL numerical flux.
 * Inputs:
 *      hmin - threshold of water depth;
 *      gra - gravity acceleration;
 *		etaM - depth of local element;
 * 		qM  - flux of local element;
 *      botM - bottom elevation of local element;
 *		etaP - depth of adjacent element;
 *		qP  - flux of adjacent element;
 *      botP - bottom elevation of adjacent element;
 *      nx  - outward normal vector;
 * Outputs:
 *      dFhs - numerical flux for water depth;
 *      dFqs - numerical flux for water flux;
 * Usages:
 * 		[dFhs, dFqs] = SWE_Mex_HLLFlux1d(hmin, gra,
 *          etaM, qM, botM, etaP, qP, botP, nx);
 */
void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    /* check input & output */
	if (nrhs != 9) mexErrMsgTxt("The number of input should be 9.");
	if (nlhs != 2) mexErrMsgTxt("The number of output should be 2");

    /* get inputs */
	double hmin = mxGetScalar(prhs[0]);
	double gra  = mxGetScalar(prhs[1]);
	double *etaM = mxGetPr(prhs[2]);
	double *qM   = mxGetPr(prhs[3]);
    double *botM = mxGetPr(prhs[4]);
    double *etaP = mxGetPr(prhs[5]);
	double *qP   = mxGetPr(prhs[6]);
    double *botP = mxGetPr(prhs[7]);
    double *nx  = mxGetPr(prhs[8]);

	size_t Nfp = mxGetM(prhs[2]);
	size_t K = mxGetN(prhs[2]);
	/* allocate output array */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);

    double *dFhs = mxGetPr(plhs[0]);
    double *dFqs = mxGetPr(plhs[1]);

    int n,k,sk=0;
    for(k=0;k<K;k++){
        for(n=0;n<Nfp;n++){
            double eM_ = etaM[sk];
            double qM_ = qM[sk];
            double bM_ = botM[sk];
            double eP_ = etaP[sk];
            double qP_ = qP[sk];
            double bP_ = botP[sk];
            double nx_ = nx[sk];
            double qnM =  qM_*nx_;
            double qnP =  qP_*nx_;

            double Fhn, Fqn;
            hll_nodal_flux(hmin, gra, eM_, eP_, qnM, qnP, bM_, bP_, &Fhn, &Fqn);
            dFhs[sk] = -Fhn;
            dFqs[sk] = -Fqn*nx_;

            nodal_flux(hmin, gra, eM_, qM_, bM_, &Fhn, &Fqn);
            dFhs[sk] += Fhn*nx_; // FhM*nx - Fhs
            dFqs[sk] += Fqn*nx_; // Fhq*nx - Fqs
            sk++;
        }
    }
    return;
}
