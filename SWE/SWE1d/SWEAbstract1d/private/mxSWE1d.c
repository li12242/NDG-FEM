#include "mxSWE1d.h"

void evaluateFlowRateByDeptheThreshold(double hcrit,
                                       double h,
                                       double hu,
                                       double* u) {
  if (h > hcrit) {
    *u = hu / h;
  } else {
    double h4 = h * h * h * h;
    *u = sqrt(2.0) * h * hu / sqrt(h4 + max(h4, TAU));
  }
  return;
}