#include "SWE2d.h"

/* Computation of the numerical flux for boundary faces.
 * Usages:
 * 	  [Fhs, Fqxs, Fqys] = SWE_Flux2d(hmin, gra, hM, hP, qxM, qxP, qyM, qyP, nx, ny)
 */
void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){

	/* check input & output */
	if (nrhs != 10)
		mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 3)
		mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	real hmin = mxGetScalar(prhs[0]);
	real gra  = mxGetScalar(prhs[1]);
	real *hM  = mxGetPr(prhs[2]);
	real *hP  = mxGetPr(prhs[3]);
	real *qxM = mxGetPr(prhs[4]);
	real *qxP = mxGetPr(prhs[5]);
	real *qyM = mxGetPr(prhs[6]);
	real *qyP = mxGetPr(prhs[7]);
	real *nx  = mxGetPr(prhs[8]);
    real *ny  = mxGetPr(prhs[9]);

	/* get dimensions */
    size_t Nfp, K;
    Nfp = mxGetM(prhs[2]);
    K   = mxGetN(prhs[2]);

	/* allocate output array */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);
	plhs[2] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);

	real *Fhs  = mxGetPr(plhs[0]); 
    real *Fqxs = mxGetPr(plhs[1]); 
    real *Fqys = mxGetPr(plhs[2]);

    int i,k,sk=0;
    for(k=0;k<K;k++){
    	for(i=0;i<Nfp;i++){
    		real hMf, hPf;
    		real qxMf, qxPf, qyMf, qyPf;
    		real qnM, qnP, qvM, qvP;
    		real nxf, nyf;
    		real Fhns, Fqns, Fqyns;

    		nxf = nx[sk];
            nyf = ny[sk];

    		hMf = hM[sk];
    		hPf = hP[sk];
    		
            qxMf = qxM[sk]; qyMf = qyM[sk];
            qxPf = qxP[sk]; qyPf = qyP[sk];

            qnM =  qxMf*nxf + qyMf*nyf;
            qvM = -qxMf*nyf + qyMf*nxf;

            qnP =  qxPf*nxf + qyPf*nyf;
            qvP = -qxPf*nyf + qyPf*nxf;

            // mexPrintf("t=%d,hM=%f,hP=%f,qnM=%f,qnP=%f,qvM=%f,qvP=%f\n",
            // 	sk,hMf,hPf,qnM,qnP,qvM,qvP);

            SWE_HLL2d(hmin, gra, hMf, hPf, qnM, qnP, qvM, qvP, 
                &Fhns, &Fqns, &Fqyns);

            Fhs [sk] = Fhns;
            Fqxs[sk] = Fqns*nxf - Fqyns*nyf;
            Fqys[sk] = Fqns*nyf + Fqyns*nxf;

            sk++;
    	}
    }

	return;
}