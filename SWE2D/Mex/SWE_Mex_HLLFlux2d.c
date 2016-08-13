#include "mex.h"
#include <math.h>

#define real double

void SWE_HLL2d(real hmin, real gra, real hM, real hP,
               real qnM, real qnP, real qvM, real qvP,
               real *Fhn, real *Fqxn, real *Fqyn);
void SWE_NodalFlux2d(real hcrit, real gra,
                     real h, real qx, real qy,
                     real *Eh, real *Eqx, real *Eqy,
                     real *Gh, real *Gqx, real *Gqy);

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

/* SWE_Mex_HLL2d
 * Calculation of the HLL numerical flux.
 * Inputs:
 *		hM  - water depth of local element
 * 		hP  - water depth of adjacent element
 *		QnM - normal flux of the local element, Qn = Qx*nx + Qy*ny
 *		QnP - normal flux of the adjacent element
 *		QvM - tangential flux of the local element, Qv = -Qx*ny + Qy*nx
 *		QvP - tangential flux of the adjacent element
 *
 * Outputs:
 * Usages:
 * 		[Fhs, Fqxs, Fqys] = SWE_Mex_HLL2d(hmin, gra, h, Qx, Qy, nx, ny, vmapM, vmapP);
 */
void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

	/* check input & output */
	if (nrhs != 9)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 3)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	real hmin = mxGetScalar(prhs[0]);
	real gra  = mxGetScalar(prhs[1]);
	real *h   = mxGetPr(prhs[2]);
	real *qx  = mxGetPr(prhs[3]);
	real *qy  = mxGetPr(prhs[4]);
	real *nx  = mxGetPr(prhs[5]);
    real *ny  = mxGetPr(prhs[6]);
    double *vmapM = mxGetPr(prhs[7]);
    double *vmapP = mxGetPr(prhs[8]);

	/* get dimensions */
    size_t nfp, ne;
    nfp = mxGetM(prhs[5]);
    ne = mxGetN(prhs[2]);

	/* allocate output array */
	plhs[0] = mxCreateDoubleMatrix((mwSize)nfp, (mwSize)ne, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)nfp, (mwSize)ne, mxREAL);
	plhs[2] = mxCreateDoubleMatrix((mwSize)nfp, (mwSize)ne, mxREAL);

	real *Fhs  = mxGetPr(plhs[0]); 
    real *Fqxs = mxGetPr(plhs[1]); 
    real *Fqys = mxGetPr(plhs[2]);

    int i,j,ind=0;
	for (i=0;i<ne;i++){
		for(j=0;j<nfp;j++){
            real hM, hP, qxM, qxP, qyM, qyP;
            real qnM, qnP, qvM, qvP;
            real nxf, nyf;
            real Fhns, Fqns, Fqyns;
            int iM = (int)vmapM[ind]-1;
            int iP = (int)vmapP[ind]-1;

            // mexPrintf("vmapM[%d] = %d, vmapP[%d] = %d\n", ind, iM, ind, iP);
            
            hM  = h[iM];
            hP  = h[iP];
            qxM = qx[iM];
            qxP = qx[iP];
            qyM = qy[iM];
            qyP = qy[iP];
            nxf = nx[ind];
            nyf = ny[ind];

            qnM =  qxM*nxf + qyM*nyf;
            qvM = -qxM*nyf + qyM*nxf;

            qnP =  qxP*nxf + qyP*nyf;
            qvP = -qxP*nyf + qyP*nxf;

			SWE_HLL2d(hmin, gra, hM, hP, qnM, qnP, qvM, qvP, 
                &Fhns, &Fqns, &Fqyns);
			
            Fhs [ind] = Fhns;
            Fqxs[ind] = Fqns*nxf - Fqyns*nyf;
            Fqys[ind] = Fqns*nyf + Fqyns*nxf;

            ind++;
		}
	}
    
    return;
}

void SWE_HLL2d(real hmin, real gra, real hM, real hP,
               real qnM, real qnP, real qvM, real qvP,
               real *Fhn, real *Fqxn, real *Fqyn){

    real EhM,EqnM,EqvM,EhP,EqnP,EqvP;
    real GhM,GqnM,GqvM,GhP,GqnP,GqvP;

    SWE_NodalFlux2d(hmin, gra, hM, qnM, qvM, &EhM, &EqnM, &EqvM, &GhM, &GqnM, &GqvM);
    SWE_NodalFlux2d(hmin, gra, hP, qnP, qvP, &EhP, &EqnP, &EqvP, &GhP, &GqnP, &GqvP);

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
        *Fhn = EhM; *Fqxn = EqnM; *Fqyn = EqvM;
    }else if((sM<0) & (sP>0)){
        *Fhn  = (sP*EhM  - sM*EhP  + sM*sP*(hP  - hM ))/(sP - sM);
        *Fqxn = (sP*EqnM - sM*EqnP + sM*sP*(qnP - qnM))/(sP - sM);
        *Fqyn = (sP*EqvM - sM*EqvP + sM*sP*(qvP - qvM))/(sP - sM);
    }else if( (sM<0)&(sP<=0) ){
        *Fhn = EhP; *Fqxn = EqnP; *Fqyn = EqvP;
    }else if( (sM==0) & (sP==0) ){
        *Fhn = 0; *Fqxn = 0; *Fqyn = 0;
    }else{
        mexErrMsgIdAndTxt("SWE_Mex_HLL2d:",
        "The wave speed computation occurs an error.");
    }
}

void SWE_NodalFlux2d(real hcrit, real gra,
                     real h, real qx, real qy,
                     real *Eh, real *Eqx, real *Eqy,
                     real *Gh, real *Gqx, real *Gqy){

    // real hcrit = (real)solver->hcrit;
    // real gra  = (real)solver->gra;

    if(h>hcrit){
        *Eh  = qx;
        *Eqx = (real)(qx*qx/h + 0.5*gra*h*h);
        *Eqy = qx*qy/h;
        *Gh  = qy;
        *Gqx = qx*qy/h;
        *Gqy = (real)(qy*qy/h + 0.5*gra*h*h);
    }else{
        *Eh  = 0; *Eqx = 0; *Eqy = 0;
        *Gh  = 0; *Gqx = 0; *Gqy = 0;
    }
}