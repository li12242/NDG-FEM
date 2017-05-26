#include "VB_TVB.h"

void weight_vertex_value(int Nvert, int MaxKv, double *Kv,
    double *VToE, double *VToC, double *hc, double *hv)
{
    int n,k;
    for(n=0;n<Nvert;n++){ // loop for all vertex
        hv[n] = 0.0;
        int Ne = (int) Kv[n];

        double *vtoe = VToE + n*MaxKv; // nth column
        double *vtoc = VToC + n*MaxKv;
        for(k=0;k<Ne;k++){ // loop for adjacent elements
            int eind = (int) vtoe[k]-1; // to C type
            double w = vtoc[k];
            hv[n] += w*hc[eind];
        }
    }
    return;
}

void minmod(int n, double *a, double *m){
	int i;
	double f, s, t;
	/* use the first element to initialize */
	t = fabs(a[0]);
	f = sign(a[0]);
	s = f;

	for(i=1;i<n;i++){
		f = sign(a[i]);

		if (s*f<0.0){
			*m = 0.0;
			return;
		}
		t = min( t, fabs(a[i]) );
	}

	*m = t*f;
	return;
}

/*
 * modified minmod function of TVB limiter.
 */
double TVB_minmod(double a1, double a2, double dx2, double factor){
    double m;
    if(fabs(a1)<factor*dx2){
        m = a1;
    }else{
        double a[2];
        a[0] = a1; a[1] = a2;
        minmod(2, a, &m);
    }
    return m;
}

void Grad2Node(int Np, double hmean,
    double xc, double yc,
 	double *x, double *y,
    double phpx, double phpy, double *h)
{
 	int i;
 	for(i=0;i<Np;i++){
            double dx = x[i] - xc;
            double dy = y[i] - yc;
            h[i] = hmean + dx*phpx + dy*phpy;
    }
}

void MatrixSolver2(double *a, double *f, double *x){
   double det = a[0]*a[3] - a[1]*a[2];
   x[0] = ( f[0]*a[3] - f[1]*a[1])/det;
   x[1] = (-f[0]*a[2] + f[1]*a[0])/det;
   return;
}

void mexFunction (int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
    /* check input & output */
    if (nrhs != 12) mexErrMsgTxt("Wrong number of input arguments.");
    if (nlhs != 1) mexErrMsgTxt("Wrong number of output arguments.");

    /* get inputs */
    double *h = mxGetPr(prhs[0]);
    double *x = mxGetPr(prhs[1]);
    double *y = mxGetPr(prhs[2]);
    double *hmean = mxGetPr(prhs[3]);
    double *xmean = mxGetPr(prhs[4]);
    double *ymean = mxGetPr(prhs[5]);
    double *Kv = mxGetPr(prhs[6]);
    double *VToE = mxGetPr(prhs[7]);
    double *VToC = mxGetPr(prhs[8]);
    double *EToV = mxGetPr(prhs[9]);
    double *Fmask = mxGetPr(prhs[10]);
    double factor = mxGetScalar(prhs[11]);

    /* get dimensions */
    size_t MaxKv = mxGetM(prhs[7]);
    size_t Nvert = mxGetN(prhs[7]);
    size_t Np = mxGetM(prhs[0]);
    size_t K = mxGetN(prhs[0]);
    size_t Nfp = mxGetM(prhs[10]);
    size_t Nv = mxGetN(prhs[10]);

    /* allocation of output */
    plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
    double *hlim = mxGetPr(plhs[0]);

    /* reconstructed vertex value */
    double *hv = calloc(Nvert, sizeof(double));
    weight_vertex_value(Nvert, MaxKv, Kv, VToE, VToC, hmean, hv);

    int k,n;
    for(k=0;k<K;k++){
        double *etov = EToV + Nv*k;
        double hc = hmean[k];
        double xc = xmean[k];
        double yc = ymean[k];
        double *h_ = h+k*Np;
        double *x_ = x+k*Np;
        double *y_ = y+k*Np;

        double delta[Nv], vx_[Nv], vy_[Nv];
        for(n=0;n<Nv;n++){
            int locid = (int)Fmask[n*Nfp] - 1; // change to C type
            double hv1 = h_[locid] - hc;

            int vind = (int)etov[n] - 1; // change to C type
            vx_[n] = x_[locid];
            vy_[n] = y_[locid];
            double dx2 = (xc-vx_[n])*(xc-vx_[n]) + (yc-vy_[n])*(yc-vy_[n]);

            double hv2 = hv[vind] - hc; // weighted vertex value
            delta[n] = TVB_minmod(hv1, 1.5*hv2, dx2, factor);
        }

        /* correct the deviations */
        double sum = 0.0;
        for(n=0;n<Nv;n++){ sum += delta[n]; }

        if( abs(sum)>TOTALERR ){
            double pos = 0.0;
            double neg = 0.0;
            for(n=0;n<Nv;n++){
                pos += max( delta[n], 0.0);
                neg += max(-delta[n], 0.0);
            }
            for(n=0;n<Nv;n++){
                delta[n] = min(1.0, neg/pos)*max(delta[n], 0.0)
                    - min(1.0, pos/neg)*max(-delta[n], 0.0);
            }
        }


        double a[4], b[2], alpha[2];
        a[0] = vx_[0] - vx_[1]; a[1] = vy_[0] - vy_[1]; b[0] = delta[0] - delta[1];
        a[2] = vx_[0] - vx_[2]; a[3] = vy_[0] - vy_[2]; b[1] = delta[0] - delta[2];
        MatrixSolver2(a, b, alpha);
        double qpx = alpha[0];
        double qpy = alpha[1];
        Grad2Node(Np, hc, xc, yc, x_, y_, qpx, qpy, hlim+k*Np);
    }
    free(hv);
    return;
}
