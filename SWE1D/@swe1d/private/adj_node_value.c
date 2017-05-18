#include "swe.h"

/*
 * Return the adjacent element boundary node values
 * according to the boundary types.
 *
 * Usages:
 *  [hP, qP] = adj_node_value(hM, qM, hP, qP, eidtype);
 */
void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    /* check input & output */
	if (nrhs != 7) mexErrMsgTxt("The number of input should be 7.");
	if (nlhs != 2) mexErrMsgTxt("The number of output should be 2.");

    /* get inputs */
	double *hM = mxGetPr(prhs[0]);
	double *qM = mxGetPr(prhs[1]);
	double *hP = mxGetPr(prhs[2]);
	double *qP = mxGetPr(prhs[3]);
    double *hext = mxGetPr(prhs[4]);
    double *qext = mxGetPr(prhs[5]);
    signed char *eidtype = (signed char *)mxGetData(prhs[6]);

    size_t Nfp = mxGetM(prhs[0]);
	size_t K = mxGetN(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);
	plhs[1] = mxCreateDoubleMatrix((mwSize)Nfp, (mwSize)K, mxREAL);

    double *hP_ = mxGetPr(plhs[0]);
    double *qP_ = mxGetPr(plhs[1]);

    int k,n,sk=0;
    for(k=0;k<K;k++){
        for(n=0;n<Nfp;n++){
            bc_type facetype = (bc_type) eidtype[sk];
            switch (facetype){
                case INNER:
                    hP_[sk] = hP[sk]; qP_[sk] = qP[sk]; break;
                case SLIPWALL:
                    hP_[sk] = hM[sk]; qP_[sk] = -qM[sk]; break;
                case NSLIPWALL:
                    hP_[sk] = hM[sk]; qP_[sk] = -qM[sk]; break;
                case ZEROGRAD:
                    hP_[sk] = hM[sk]; qP_[sk] = qM[sk]; break;
                case CLAMP:
                    hP_[sk] = hext[sk]; qP_[sk] = qext[sk]; break;
            }
            sk++;
        }
    }
    return;
}
