#include "swe.h"

void hll_nodal_flux(double hmin, double gra, double hM, double hP,
    double qnM, double qnP, double *Fhn, double *Fqxn){

    double FhM, FqnM, FhP, FqnP;
    reduce_nodal_flux(hmin, gra, hM, qnM, &FhM, &FqnM);
    reduce_nodal_flux(hmin, gra, hP, qnP, &FhP, &FqnP);

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
        mexPrintf("[hM, qnM]=[%f,%f], [hP, qnP]=[%f,%f], [sM,sP]=[%f,%f]",
                hM, qnM, hP, qnP, sM, sP);
        mexErrMsgIdAndTxt("MATLAB:hll_flux:hll_nodal_flux:",
        "The wave speed computation occurs an error.");
    }
    return;
}

/* SWE_Mex_HLL2d
 * Calculation of the HLL numerical flux.
 * Inputs:
 *      hmin - threshold of water depth;
 *      gra - gravity acceleration;
 *		hM  - depth of local element;
 * 		qM  - flux of local element;
 *		hP  - depth of adjacent element;
 *		qP  - flux of adjacent element;
 *      nx  - outward normal vector;
 * Outputs:
 *      dFhs - numerical flux for water depth;
 *      dFqs - numerical flux for water flux;
 * Usages:
 * 		[dFhs, dFqs] = SWE_Mex_HLLFlux1d(hmin, gra, hM, qM, hP, qP, zM, nx);
 */
void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    /* check input & output */
	if (nrhs != 7) mexErrMsgTxt("The number of input should be 7.");
	if (nlhs != 2) mexErrMsgTxt("The number of output should be 2");

    /* get inputs */
	double hmin = mxGetScalar(prhs[0]);
	double gra  = mxGetScalar(prhs[1]);
	double *hM  = mxGetPr(prhs[2]);
	double *qM  = mxGetPr(prhs[3]);
    double *hP  = mxGetPr(prhs[4]);
	double *qP  = mxGetPr(prhs[5]);
    double *nx  = mxGetPr(prhs[6]);

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
            double hM_ = hM[sk];
            double qM_ = qM[sk];
            double hP_ = hP[sk];
            double qP_ = qP[sk];
            double nx_ = nx[sk];
            double qnM =  qM_*nx_;
            double qnP =  qP_*nx_;

            double Fhn, Fqn;
            hll_nodal_flux(hmin, gra, hM_, hP_, qnM, qnP, &Fhn, &Fqn);
            dFhs[sk] = -Fhn;
            dFqs[sk] = -Fqn*nx_;

            reduce_nodal_flux(hmin, gra, hM_, qM_, &Fhn, &Fqn);
            dFhs[sk] += Fhn*nx_; // FhM*nx - Fhs
            dFqs[sk] += Fqn*nx_; // Fhq*nx - Fqs
            sk++;
        }
    }
    return;
}
