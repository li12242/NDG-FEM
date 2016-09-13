#include "SWE2d.h"

/**
 * @brief Calculate the flux term of each node
 * @param [real] hcrit  threadhold of the water depth
 * @param [real] gra    gravity accelerated
 * @param [real] h      water depth
 * @param [real] qx     water flux along x coordinate
 * @param [real] qy     water flux along y coordinate
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * E | real*   | the flux terms along x coordinate
 * G | real*   | the flux terms along y coordinate
 *
 */
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

/**
 * @brief Calculation of the numerical flux of each node
 * @param [real] hmin   threadhold of the water depth
 * @param [real] gra    gravity accelerated
 * @param [real] hM     water depth on local node
 * @param [real] hP     water depth on adjacent node
 * @param [real] qnM    
 * @param [real] qnP
 * @param [real] qvM
 * @param [real] qvP
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * Fhn  | real*   | the flux term of variable h
 * Fqxn | real*   | the flux term of variable qx
 * Fqyn | real*   | the flux term of variable qy
 * 
 */
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
