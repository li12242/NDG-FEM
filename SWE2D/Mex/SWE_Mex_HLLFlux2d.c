#include "SWE2d.h"

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
    real *vmapM = mxGetPr(prhs[7]);
    real *vmapP = mxGetPr(prhs[8]);

	/* get dimensions */
    size_t nfp, ne;
    nfp = mxGetM(prhs[5]);
    ne  = mxGetN(prhs[2]);

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

