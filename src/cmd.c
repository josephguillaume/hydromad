// hydromad: Hydrological Modelling and Analysis of Data
//
// Copyright (c) Felix Andrews <felix@nfrac.org>

#include <R_ext/Arith.h>

#include <R.h>

#ifndef min
#define min(a, b) ((a < b) ? (a) : (b))
#define max(a, b) ((a < b) ? (b) : (a))
#endif

void sma_cmd(double *P, double *E, int *n, double *g, double *e, double *d,
             double *shape, double *M_0, double *U, double *M, double *ET) {
  int t;
  double Mf, a;
  double M_prev = *M_0;
  for (t = 0; t < *n; t++) {
    // default, for when P[t] == 0:
    Mf = M_prev;
    // select form of dU/dP relationship
    if (P[t] > 0) {
      // rainfall reduces CMD (Mf)
      if (*shape < 1) {
        // linear form: dU/dP = 1 - (M/d)
        if (M_prev < *d) {
          Mf = M_prev * exp(-P[t] / *d);
        } else if (M_prev < *d + P[t]) {
          Mf = *d * exp((-P[t] + M_prev - *d) / *d);
        } else {
          Mf = M_prev - P[t];
        }
      } else if (*shape == 1) {
        // hyperbolic form: dU/dP = 1 - ??
        if (M_prev < *d) {
          Mf = 1 / tan((M_prev / *d) * (PI / 2));
          Mf = (2 * *d / PI) * atan(1 / (PI * P[t] / (2 * *d) + Mf));
        } else if (M_prev < *d + P[t]) {
          Mf = (2 * *d / PI) * atan(2 * *d / (PI * (*d - M_prev + P[t])));
        } else {
          Mf = M_prev - P[t];
        }
      } else { // *shape > 1
        // power form: dU/dP = 1 - (M/d)^b
        a = pow(10, (*shape / 50));
        if (M_prev < *d) {
          Mf = M_prev *
               pow(1 - ((1 - a) * P[t] / pow(*d, a)) / pow(M_prev, (1 - a)),
                   1 / (1 - a));
        } else if (M_prev < *d + P[t]) {
          Mf = *d * pow(1 - (1 - a) * (P[t] - M_prev + *d) / *d, 1 / (1 - a));
        } else {
          Mf = M_prev - P[t];
        }
      }
    }
    // drainage (rainfall not accounted for in -dM)
    U[t] = max(0, P[t] - M_prev + Mf);
    // evapo-transpiration
    ET[t] = *e * E[t] * min(1, exp(2 * (1 - Mf / *g)));
    ET[t] = max(0, ET[t]);
    // mass balance
    M[t] = M_prev - P[t] + U[t] + ET[t];
    M[t] = max(0, M[t]);
    M_prev = M[t];
  }
}
