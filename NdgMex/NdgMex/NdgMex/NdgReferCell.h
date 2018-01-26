//
//  NdgReferElement.h
//  NdgMex
//
//  Created by li12242 on 18/1/20.
//  Copyright (c) 2018年 li12242. All rights reserved.
//

#ifndef __NdgMex__NdgReferElement__
#define __NdgMex__NdgReferElement__

#include "NdgMath.h"
#include <stdio.h>

/** Reference element type */
typedef enum {
  NdgPoint = 0,
  NdgLine = 1,
  NdgTri = 2,
  NdgQuad = 3,
  NdgPrismTri = 4,
  NdgPrismQuad = 5,
} NdgTypeReferCell;

/** Reference element class */
class NdgReferCell {

private:
  int N;
  int Np;
  int Nq;
  int Nv;
  int Nface;
  int *Nfv;
  int *Nfp;
  int totalNfp;
  NdgTypeReferCell type;
  int **Fmask;
  /** coordinates of vertex, quadrature and interpolation nodes */
  NdgMath::NdgDouble *vr, *vs, *vt;
  NdgMath::NdgDouble *r, *s, *t;
  NdgMath::NdgDouble *rq, *sq, *tq, *wq;
  NdgMath::NdgDouble LAV;
  int **FToV;
  NdgMath::NdgDouble **Dr, **Ds, **Dt;
  NdgMath::NdgDouble **M, **invM;
  NdgMath::NdgDouble **V, **invV;

public:
  inline int getN() { return N; };
  inline int getNp() { return Np; };
  inline int getNv() { return Nv; };
  inline int getNq() { return Nq; };
  inline int getTotalNfp() { return totalNfp; };
  inline int getNfp(int n) { return Nfp[n]; };
  inline NdgTypeReferCell getCellType() { return type; };
};

#endif /* defined(__NdgMex__NdgReferElement__) */
