#include "mex.h"
#include <math.h>

#define real double

void SWE_HLL1d(real hmin, real gra, real hM, real hP,
               real qnM, real qnP, real *Fhn, real *Fqxn);

void SWE_NodalFlux1d(real hcrit, real gra,
                     real h,   real qx,
                     real *Fh, real *Fq);

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

/* SWE_Mex_HLL2d
 * Calculation of the HLL numerical flux.
 * Inputs:
 *		hM  - water depth of local element
 * 		hP  - water depth of adjacent element
 *		QnM - normal flux of the local element, Qn = Qx*nx + Qy*ny
 *		QnP - normal flux of the adjacent element
 *
 * Outputs:
 * Usages:
 * 		[Fhs, Fqxs] = SWE_Mex_HLLFlux1d(hmin, gra, h, Qx, nx, vmapM, vmapP);
 */
void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

	/* check input & output */
	if (nrhs != 7)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 2)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	real hmin = mxGetScalar(prhs[0]);
	real gra  = mxGetScalar(prhs[1]);
	real *h   = mxGetPr(prhs[2]);
	real *qx  = mxGetPr(prhs[3]);
	real *nx  = mxGetPr(prhs[4]);
    double *vmapM = mxGetPr(prhs[5]);
    double *vmapP = mxGetPr(prhs[6]);

	/* get dimensions */
    size_t nfp, ne;
    nfp = mxGetM(prhs[5]);
    ne  = mxGetN(prhs[2]);

	/* allocate output array */
	plhs[0] = mxCreateDoubleMatrix((mwSize)nfp, (mwSize)ne, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)nfp, (mwSize)ne, mxREAL);

	real *Fhs  = mxGetPr(plhs[0]); 
    real *Fqxs = mxGetPr(plhs[1]); 

    int i,j,ind=0;
	for (i=0;i<ne;i++){
		for(j=0;j<nfp;j++){
            real hM, hP, qxM, qxP;
            real qnM, qnP;
            real nxf;
            real Fhns, Fqns;

            int iM = (int)vmapM[ind]-1;
            int iP = (int)vmapP[ind]-1;

            // mexPrintf("vmapM[%d] = %d, vmapP[%d] = %d\n", ind, iM, ind, iP);
            
            hM  = h[iM];
            hP  = h[iP];
            qxM = qx[iM];
            qxP = qx[iP];
            nxf = nx[ind];

            qnM =  qxM*nxf;
            qnP =  qxP*nxf;

			SWE_HLL1d(hmin, gra, hM, hP, qnM, qnP,
                &Fhns, &Fqns);
			
            Fhs [ind] = Fhns;
            Fqxs[ind] = Fqns*nxf;
            ind++;
		}
	}
    
    return;
}


void SWE_HLL1d(real hmin, real gra, real hM, real hP,
               real qnM, real qnP, real *Fhn, real *Fqxn){

    real FhM,FqnM,FhP,FqnP;

    SWE_NodalFlux1d(hmin, gra, hM, qnM, &FhM, &FqnM);
    SWE_NodalFlux1d(hmin, gra, hP, qnP, &FhP, &FqnP);

    /* calculation of wave speed */
    real sM, sP;
    // real gra  = (real) solver->gra;
    // real hmin = (real) solver->hcrit;
    real us, cs;
    real unM,unP;

    if( (hM>hmin) & (hP>hmin) ){
        unM=qnM/hM;
        unP=qnP/hP;
        us = (real)(0.5*(unM + unP)   + sqrt(gra*hM)   - sqrt(gra*hP));
        cs = (real)(0.5*(sqrt(gra*hM) + sqrt(gra*hP) ) + 0.25*(unM - unP));

        sM = (real)min(unM-sqrt(gra*hM), us-cs);
        sP = (real)max(unP+sqrt(gra*hP), us+cs);
    }else if ( (hM>hmin) & (hP<=hmin) ){
        unM=qnM/hM;
        sM = (real)(unM -  sqrt(gra*hM) );
        sP = (real)(unM +2*sqrt(gra*hM) );
    }else if ( (hM<=hmin) & (hP>hmin) ){
        unP=qnP/hP;
        sM = (real)(unP -2*sqrt(gra*hP) );
        sP = (real)(unP +  sqrt(gra*hP) );
    }else{ /* both dry element */
        sM = 0; sP = 0;
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

void SWE_NodalFlux1d(real hcrit, real gra,
                     real h,   real qx,
                     real *Fh, real *Fq){
    if(h>hcrit){
        *Fh  = qx;
        *Fq = (real)(qx*qx/h + 0.5*gra*h*h);
    }else{
        *Fh  = 0; *Fq = 0;
    }

    return;
}