// hydromad: Hydrological Modelling and Analysis of Data
//
// Copyright (c) Felix Andrews <felix@nfrac.org>

#include <R_ext/Arith.h>

#include <R.h>

#ifndef min
#define min(a, b) ((a < b) ? (a) : (b))
#define max(a, b) ((a < b) ? (b) : (a))
#endif

#define my_isok(x) (!ISNA(x) & !ISNAN(x))

void ar1_tv(double *x, double *a, int *n, double *init, double *out) {
  int i;
  double tmp;
  tmp = *init;
  for (i = 0; i < *n; i++) {
    out[i] = x[i] + a[i] * tmp;
    tmp = out[i];
  }
}

void filter_constloss(double *x, int *n, double *a, int *na, double *loss,
                      double *out) {
  int i, j;
  // (skip over first *na elements)
  for (i = *na; i < *n; i++) {
    out[i] = x[i];
    for (j = 0; j < *na; j++) {
      out[i] = out[i] + a[j] * out[i - j - 1];
    }
    out[i] = max(0, out[i] - *loss);
  }
}
