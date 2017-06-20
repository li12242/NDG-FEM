#include "vb.h"
#define EPSILON 1.0e-12

/* the weights of van Albada limiter */
void VAGrad(int Nsub, double *gra_x, double *gra_y, double *gra_det,
    double *dhdx, double *dhdy){
    double frac=Nsub*EPSILON;;
    int i,j;

    for(*dhdx=0.0,*dhdy=0.0,i=0;i<Nsub;i++){
        double w = 1.0;
        for(j=0;j<Nsub;j++){
            if(i==j)
                continue;
            w = w*gra_det[j];
        }
        w += EPSILON;
        frac += w;
        *dhdx += w*gra_x[i];
        *dhdy += w*gra_y[i];
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
        Fmask, EToV, flim, VAGrad);

    return;
}
