#include "swe.h"
#define DEBUG 0
/**
 * @brief Calculation of the Lax-Friedrichs numerical flux
 * @param [in] hmin - threadhold of the water depth
 * @param [in] gra - gravity accelerated
 * @param [in] hM, qnM, qvM - water depth on local element node
 * @param [in] hP, qnP, qvP - water depth on adjacent element node
 * @param [out] Fhn, Fqxn, Fqyn - HLL numerical flux;
 */
void lf_flux(double hmin, double gra, double hM, double hP,
              double qnM, double qnP, double qvM, double qvP,
              double *Fhn, double *Fqxn, double *Fqyn)
{
    double EhM,EqnM,EqvM,EhP,EqnP,EqvP;
    double GhM,GqnM,GqvM,GhP,GqnP,GqvP;
    nodal_flux(hmin, gra, hM, qnM, qvM, &EhM, &EqnM, &EqvM, &GhM, &GqnM, &GqvM);
    nodal_flux(hmin, gra, hP, qnP, qvP, &EhP, &EqnP, &EqvP, &GhP, &GqnP, &GqvP);

    double sM, sP;
    if( (hM>hmin) ){ sM = fabs(qnM/hM) + sqrt(gra*hM);
    }else{ sM = 0.0; }
    if( (hP>hmin) ){ sP = fabs(qnP/hP) + sqrt(gra*hP);
    }else{ sP = 0.0; }

    #if DEBUG
    mexPrintf("hM=%f, hP=%f, qnM=%f, qnP=%f, sM=%f, sP=%f\n",
        hM, hP, qnM, qnP, sM, sP);
    #endif

    double s = max(sM, sP);
    *Fhn = 0.5*( EhM + EhP + s*(hM-hP) );
    *Fqxn = 0.5*( EqnM + EqnP + s*(qnM-qnP) );
    *Fqyn = 0.5*( EqvM + EqvP + s*(qvM-qvP) );

    #if DEBUG
    mexPrintf("Fhn=%f, Fqxn=%f, Fqyn=%f\n", *Fhn, *Fqxn, *Fqyn);
    #endif
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
 *      h, qx, qy, h_ext, qx_ext, qy_ext, nx, ny, eidM, eidP, eidtype);
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
    signed char *eidtype = (signed char *)mxGetData(prhs[12]); // int8 类型

	/* get dimensions */
    size_t Nfp = mxGetM(prhs[10]);
    size_t K = mxGetN(prhs[10]);

	/* allocate output array */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);
	plhs[2] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);

	double *dFh  = mxGetPr(plhs[0]);
    double *dFqx = mxGetPr(plhs[1]);
    double *dFqy = mxGetPr(plhs[2]);

    int i,j,ind=0;
	for (i=0;i<K;i++){
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
            if(info) mexErrMsgTxt("Unknown boundary conditions.");

            double qnM, qnP, qvM, qvP;
            qnM = f_M[1]*nx_ + f_M[2]*ny_; qvM = -f_M[1]*ny_ + f_M[2]*nx_;
            qnP = f_P[1]*nx_ + f_P[2]*ny_; qvP = -f_P[1]*ny_ + f_P[2]*nx_;
            //qnM =  qxM*nx_ + qyM*ny_; qvM = -qxM*ny_ + qyM*nx_;
            //qnP =  qxP*nx_ + qyP*ny_; qvP = -qxP*ny_ + qyP*nx_;

            double Fhns, Fqns, Fqyns;
		    lf_flux(hmin, gra, f_M[0], f_P[0],
                qnM, qnP, qvM, qvP, &Fhns, &Fqns, &Fqyns);

            dFh[ind] = -Fhns;
            dFqx[ind] = -(Fqns*nx_ - Fqyns*ny_);
            dFqy[ind] = -(Fqns*ny_ + Fqyns*nx_);

            double Eh, Eqx, Eqy, Gh, Gqx, Gqy;
            nodal_flux(hmin, gra, f_M[0], f_M[1], f_M[2],
                &Eh, &Eqx, &Eqy, &Gh, &Gqx, &Gqy);

            dFh[ind] += nx_*Eh + ny_*Gh;
            dFqx[ind] += nx_*Eqx + ny_*Gqx;
            dFqy[ind] += nx_*Eqy + ny_*Gqy;

            ind++;
		}
	}

    return;
}
