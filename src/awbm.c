// hydromad: Hydrological Modelling and Analysis of Data
//
// Copyright (c) Felix Andrews <felix@nfrac.org>

#include <R_ext/Arith.h>

#include <R.h>

#ifndef min
#define min(a, b) ((a < b) ? (a) : (b))
#define max(a, b) ((a < b) ? (b) : (a))
#endif

void sma_awbm(double *P, double *E, int *n, double *cap1, double *cap2,
              double *cap3, double *area1, double *area2, double *area3,
              double *S1_0, double *S2_0, double *S3_0, double *U, double *S1,
              double *S2, double *S3) {
  int t;
  double U1, U2, U3;
  double S1_prev = *S1_0;
  double S2_prev = *S2_0;
  double S3_prev = *S3_0;
  for (t = 0; t < *n; t++) {
    // rainfall
    S1[t] = S1_prev + P[t];
    S2[t] = S2_prev + P[t];
    S3[t] = S3_prev + P[t];
    // saturation excess
    U1 = max(0, S1[t] - *cap1);
    U2 = max(0, S2[t] - *cap2);
    U3 = max(0, S3[t] - *cap3);
    U[t] = *area1 * U1 + *area2 * U2 + *area3 * U3;
    // evapotranspiration
    S1_prev = S1[t] = max(0, S1[t] - U1 - E[t]);
    S2_prev = S2[t] = max(0, S2[t] - U2 - E[t]);
    S3_prev = S3[t] = max(0, S3[t] - U3 - E[t]);
  }
}
