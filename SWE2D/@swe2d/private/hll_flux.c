#include "swe.h"

/**
 * @brief Calculation of the HLL numerical flux
 * @param [in] hmin - threadhold of the water depth
 * @param [in] gra - gravity accelerated
 * @param [in] hM, qnM, qvM - water depth on local element node
 * @param [in] hP, qnP, qvP - water depth on adjacent element node
 * @param [out] Fhn, Fqxn, Fqyn - HLL numerical flux;
 */
void hll_flux(double hmin, double gra, double hM, double hP,
              double qnM, double qnP, double qvM, double qvP,
              double *Fhn, double *Fqxn, double *Fqyn){

    double EhM,EqnM,EqvM,EhP,EqnP,EqvP;
    double GhM,GqnM,GqvM,GhP,GqnP,GqvP;
    reduce_nodal_flux(hmin, gra, hM, qnM, qvM,
        &EhM, &EqnM, &EqvM, &GhM, &GqnM, &GqvM);
    reduce_nodal_flux(hmin, gra, hP, qnP, qvP,
        &EhP, &EqnP, &EqvP, &GhP, &GqnP, &GqvP);
    /* calculation of wave speed */
    double sM, sP, us, cs, unM,unP;
    if( (hM>hmin) & (hP>hmin) ){
        unM=qnM/hM;
        unP=qnP/hP;
        us = (double)(0.5*(unM + unP)   + sqrt(gra*hM)   - sqrt(gra*hP));
        cs = (double)(0.5*(sqrt(gra*hM) + sqrt(gra*hP) ) + 0.25*(unM - unP));

        sM = (double)min(unM-sqrt(gra*hM), us-cs);
        sP = (double)max(unP+sqrt(gra*hP), us+cs);
    }else if ( (hM>hmin) & (hP<=hmin) ){
        unM=qnM/hM;
        sM = (double)(unM -  sqrt(gra*hM) );
        sP = (double)(unM +2*sqrt(gra*hM) );
    }else if ( (hM<=hmin) & (hP>hmin) ){
        unP=qnP/hP;
        sM = (double)(unP -2*sqrt(gra*hP) );
        sP = (double)(unP +  sqrt(gra*hP) );
    }else{ /* both dry element */
        sM = 0; sP = 0;
    }
    /* HLL flux function */
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
    return;
}

/*
 * Calculation of the flux deviations with HLL numerical flux.
 * Inputs:
 *  hmin        - minimum water depth;
 *  gra         - gravity acceleration;
 *  h, qx, qy   - conservative variables;
 *	h_ext, qx_ext, qy_ext - external values of each nodes;
 *  nx, ny      - outward normal vectors of local element;
 *  eidM, eidP  - node numbers of local and adjacent elements;
 *  eidtype     - boundary types with int8 type variable;
 *
 * Usages:
 * 	[dFhs, dFqxs, dFqys] = SWE_Mex_HLL2d(hmin, gra,
 *      h, qx, qy, z, h_ext, qx_ext, qy_ext, nx, ny, eidM, eidP, eidtype);
 */
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]){

	/* check input & output */
	if (nrhs != 13) mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 3) mexErrMsgTxt("Wrong number of output arguments.");

	/* get inputs */
	double hmin = mxGetScalar(prhs[0]);
	double gra = mxGetScalar(prhs[1]);
	double *h = mxGetPr(prhs[2]);
	double *qx = mxGetPr(prhs[3]);
	double *qy = mxGetPr(prhs[4]);
    double *h_ext = mxGetPr(prhs[5]);
	double *qx_ext = mxGetPr(prhs[6]);
	double *qy_ext = mxGetPr(prhs[7]);
	double *nx = mxGetPr(prhs[8]);
    double *ny = mxGetPr(prhs[9]);
    double *eidM = mxGetPr(prhs[10]);
    double *eidP = mxGetPr(prhs[11]);
    signed char *eidtype = (signed char *)mxGetData(prhs[12]); // int8 ç±»å?

	/* get dimensions */
    size_t Nfp = mxGetM(prhs[11]);
    size_t K = mxGetN(prhs[11]);

	/* allocate output array */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);
	plhs[2] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);

	double *dFh  = mxGetPr(plhs[0]);
    double *dFqx = mxGetPr(plhs[1]);
    double *dFqy = mxGetPr(plhs[2]);

    /* set number of threads */
//     int n = omp_get_num_procs();
//     omp_set_num_threads(n);
    int i,j;
    //#pragma omp parallel for private(j)
	for (i=0;i<K;i++){
        int ind = i*Nfp;
		for(j=0;j<Nfp;j++){
            int iM = (int)eidM[ind]-1; // change to C type
            int iP = (int)eidP[ind]-1;
            double f_M[3], varP[3]; // local and adjacent node values
            f_M[0] = h[iM];  varP[0] = h[iP];
            f_M[1] = qx[iM]; varP[1] = qx[iP];
            f_M[2] = qy[iM]; varP[2] = qy[iP];

            // outward normal vector of local element
            double nx_ = nx[ind];
            double ny_ = ny[ind];

            double f_ext[3]; // external values on local nodes
            f_ext[0] = h_ext[iM];
            f_ext[1] = qx_ext[iM];
            f_ext[2] = qy_ext[iM];

            bc_type type = (bc_type)eidtype[ind];
            // get adjacent values hP, qxP, qyP, considering
            // various boudnary conditions
            double f_P[3];
            int info = bound_cond(f_M, varP, f_ext, nx_, ny_, type, f_P);
            // if(info) mexErrMsgTxt("Unknown boundary conditions.");

            double qnM, qnP, qvM, qvP;
            qnM = f_M[1]*nx_ + f_M[2]*ny_; qvM = -f_M[1]*ny_ + f_M[2]*nx_;
            qnP = f_P[1]*nx_ + f_P[2]*ny_; qvP = -f_P[1]*ny_ + f_P[2]*nx_;
            //qnM =  qxM*nx_ + qyM*ny_; qvM = -qxM*ny_ + qyM*nx_;
            //qnP =  qxP*nx_ + qyP*ny_; qvP = -qxP*ny_ + qyP*nx_;

            double Fhns, Fqns, Fqyns;
			hll_flux(hmin, gra, f_M[0], f_P[0],
                qnM, qnP, qvM, qvP, &Fhns, &Fqns, &Fqyns);

            dFh[ind] = -Fhns;
            dFqx[ind] = -(Fqns*nx_ - Fqyns*ny_);
            dFqy[ind] = -(Fqns*ny_ + Fqyns*nx_);

            double Eh, Eqx, Eqy, Gh, Gqx, Gqy;
            reduce_nodal_flux(hmin, gra, f_M[0], f_M[1], f_M[2],
                &Eh, &Eqx, &Eqy, &Gh, &Gqx, &Gqy);

            dFh[ind] += nx_*Eh + ny_*Gh;
            dFqx[ind] += nx_*Eqx + ny_*Gqx;
            dFqy[ind] += nx_*Eqy + ny_*Gqy;

            ind++;
		}
	}

    return;
}
