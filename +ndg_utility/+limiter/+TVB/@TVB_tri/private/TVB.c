#include "TVB.h"
#define DEBUG 0

/**
 * @brief TVB limiter
 *
 * Usages:
 * 		hlim  = (h, hc, xc, yc, J, sJ, w, Mes, EToE, factor)
 */
void mexFunction (int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
    /* check input & output */
    if (nrhs != 11)
        mexErrMsgTxt("Wrong number of input arguments.");
    if (nlhs != 1)
        mexErrMsgTxt("Wrong number of output arguments");

    /* get inputs */
    double *h = mxGetPr(prhs[0]);
    double *x = mxGetPr(prhs[1]);
    double *y = mxGetPr(prhs[2]);
    double *hmean = mxGetPr(prhs[3]);
    double *xmean = mxGetPr(prhs[4]);
    double *ymean = mxGetPr(prhs[5]);
    double *hfm = mxGetPr(prhs[6]);
    double *xfm = mxGetPr(prhs[7]);
    double *yfm = mxGetPr(prhs[8]);
    double *EToE = mxGetPr(prhs[9]);
    double factor = mxGetScalar(prhs[10]);

    /* get dimensions */
    size_t Np = mxGetM(prhs[0]);
    size_t K = mxGetN(prhs[0]);
    size_t Nfaces = mxGetM(prhs[6]);

    if(Nfaces!=3){
        mexPrintf("Wrong number of Nfaces=%d for TVB_tri limiter", Nfaces);
        exit(-1);
    }

    /* allocation of output */
    plhs[0] = mxCreateDoubleMatrix((mwSize)Np, (mwSize)K, mxREAL);
    double *hlim = mxGetPr(plhs[0]);

    double alpha[2], a[4], b[2];
    int k,f,f1,f2,sf=0;

    for(k=0;k<K;k++){
        double xc = xmean[k];
        double yc = ymean[k];
        double hc = hmean[k];
        double xf[3], yf[3], hf[3]; // value of middle point on face
        double xe[3], ye[3], he[3]; // value of adjacent element
        double delta[3]; // limited middle point value
        #if DEBUG
        mexPrintf("k=%d, xc=%f, yc=%f, hc=%f\n", k, xc, yc, hc);
        #endif

        for(f=0;f<Nfaces;f++){
            hf[f] = hfm[sf];
            xf[f] = xfm[sf];
            yf[f] = yfm[sf];
            sf++;
            /* value of adjacent elements */
            int e1 = (int)EToE[k*Nfaces+f] - 1;
            if(e1 == k){ // for boundary
                he[f] = hmean[e1];
                xe[f] = 2*xf[f] - xc;
                ye[f] = 2*yf[f] - yc;
            }else{
                he[f] = hmean[e1];
                xe[f] = xmean[e1];
                ye[f] = ymean[e1];
            }
            #if DEBUG
            mexPrintf("k=%d, f=%d, e1=%d, xf=%f, yf=%f, hf=%f, xe=%f, ye=%f, he=%f\n", 
                        k, f, e1, xf[f], yf[f], hf[f], xe[f], ye[f], he[f]);
            #endif
        }

        for(f1=0;f1<Nfaces;f1++){
            /* find best other face */
            double max_alpha = -1.0;
            int    best_face = -1;
            double best_alpha[2];

            a[0] = xe[f1] - xc;
            a[2] = ye[f1] - yc;
            b[0] = xf[f1] - xc;
            b[1] = yf[f1] - yc;

            for(f2=0;f2<Nfaces;f2++){
                if(f1==f2){ continue; }

                /* calculate of alpha */
                a[1] = xe[f2] - xc;
                a[3] = ye[f2] - yc;

                MatrixSolver2(a, b, alpha);

                double alpha_det = sqrt(alpha[0]*alpha[0] + alpha[1]*alpha[1]);
                if ((alpha[0] > -TOTALERR) & (alpha[1] > -TOTALERR)
                    & alpha[0]/alpha_det > max_alpha ){

                    best_face = f2;
                    max_alpha = alpha[0]/alpha_det;
                    best_alpha[0] = alpha[0];
                    best_alpha[1] = alpha[1];
                }
            }
            #if DEBUG
            mexPrintf("k=%d, f=%d, f2=%d, alpha=[%f, %f]\n",
            	k,f1,best_face,best_alpha[0],best_alpha[1]);
            #endif
            double dhe, dhf, dx2;
            dhe = best_alpha[0]*(he[f1] - hc) + best_alpha[1]*(he[best_face] - hc);
            dhf = hf[f1] - hc;
            dx2 = (xe[f1]-xc)*(xe[f1]-xc)+(ye[f1]-yc)*(ye[f1]-yc);
            delta[f1] = TVB_minmod(dhf, 1.5*dhe, dx2, factor);
            #if DEBUG
            mexPrintf("k=%d, f=%d, du=%f, due=%f, M=%f, delta=%f\n",
            	k,f1,dhf,dhe,factor*dx2,delta[f1]);
            #endif
        }
        /* correct the average to sum to 0.0 */
        double sum = 0.0;
        for(f=0;f<Nfaces;f++){
            sum += delta[f];
        }

        if( abs(sum)>TOTALERR ){
            double pos = 0.0;
            double neg = 0.0;

            for(f=0;f<Nfaces;f++){
                pos += max( delta[f], 0.0);
                neg += max(-delta[f], 0.0);
            }

            for(f=0;f<Nfaces;f++){
                delta[f] = min(1.0, neg/pos)*max(delta[f], 0.0)
                    - min(1.0, pos/neg)*max(-delta[f], 0.0);
            }
            #if DEBUG
            mexPrintf("k=%d, pos=%e, neg=%e, delta=%f\n", k, pos, neg, delta);
            #endif
        }
        for(f=0;f<Nfaces;f++){ delta[f] += hmean[k]; }
        /* reconstruct the cell value */
        double qpx = 0.0;
        double qpy = 0.0;
        a[0] = xf[0] - xf[1]; a[1] = yf[0] - yf[1]; b[0] = delta[0] - delta[1];
        a[2] = xf[0] - xf[2]; a[3] = yf[0] - yf[2]; b[1] = delta[0] - delta[2];
        MatrixSolver2(a, b, alpha);
        qpx = alpha[0];
        qpy = alpha[1];
        Grad2Node(Np, hc, xc, yc, x+k*Np, y+k*Np, qpx, qpy, hlim+k*Np);
    }
    return;
}

/*
* modified minmod function of TVB limiter.
* m = TVB_minmod(a1, a2)
*/
double TVB_minmod(double a1, double a2, double dx2, double TVB_factor){
    double m;
    // mexPrintf("a1=%f, a2=%f, Th=%f\n", a1, a2, TVB_factor*dx2);
    if(fabs(a1)<TVB_factor*dx2){
        m = a1;
        // mexPrintf("a1=%d, Th=%f\n", fabs(a1), TVB_factor*dx2);
    }else{
        double a[2];
        a[0] = a1; a[1] = a2;
        minmod(2, a, &m);
        // mexPrintf("a1=%f, a2=%f, minmod=%f\n", a[0], a[1], m);
    }
    return m;
}
