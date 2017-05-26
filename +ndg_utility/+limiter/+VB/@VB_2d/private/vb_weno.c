#include "vb.h"
#define EPSILON 1.0e-12
#define DEBUG 0
/* weights of Hermite WENO limiter */
void WenoGrad(int Nsub, double *gra_x, double *gra_y, double *gra_det,
    double *dhdx, double *dhdy){
    int i;
    double frac = 0.0;
    double r =2.0; // a positive number
    *dhdx=0.0,*dhdy=0.0;
    for(i=0;i<Nsub;i++){
        double w = pow(sqrt(gra_det[i])+EPSILON, -r);
        frac += w;
        *dhdx += w*gra_x[i];
        *dhdy += w*gra_y[i];

        #if DEBUG
        mexPrintf("Nfaces=%d, w=%f\n", i, w);
        #endif
    }
    *dhdx /= frac;
    *dhdy /= frac;
}

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	/* check input & output */
	if (nrhs != 11) mexErrMsgTxt("Wrong number of input arguments.");
	if (nlhs != 1) mexErrMsgTxt("Wrong number of output arguments");

	/* get inputs */
	double *f = mxGetPr(prhs[0]);
	double *x = mxGetPr(prhs[1]);
	double *y = mxGetPr(prhs[2]);
	double *fc = mxGetPr(prhs[3]);
	double *xc = mxGetPr(prhs[4]);
	double *yc = mxGetPr(prhs[5]);
    double *fv = mxGetPr(prhs[6]);
	double *f_max = mxGetPr(prhs[7]);
	double *f_min = mxGetPr(prhs[8]);
	double *EToV = mxGetPr(prhs[9]);
	double *Fmask= mxGetPr(prhs[10]);

	/* get dimensions */
	size_t Np = mxGetM(prhs[0]);
	size_t K  = mxGetN(prhs[0]);
	size_t Nfp = mxGetM(prhs[10]);
	size_t Nfaces = mxGetN(prhs[10]);

	/* allocation of output */
	plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
	double *flim = mxGetPr(plhs[0]);

	VertLimit(K, Np, Nfaces, Nfp, fv, f_max, f_min,
        fc, xc, yc, f, x, y,
		Fmask, EToV, flim, WenoGrad);

	return;
}
