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

void sriv_system(double *U, double *Y, double *X, int *len, int *warmup,
                 int *order, int *delay, double *xz, double *xy, double *xx) {
  // z: regressors { Qf[k-1] ... Qf[k-n], Uf[k] ... Uf[k-m] }
  // x: instrumental variable { Xf[k-1] ... Xf[k-n], Uf[k] ... Uf[k-m] }
  // form the system E(xz') %*% theta == E(xy)
  // i.e. (xz) compute mean of outer products of x_t and z_t for all times t
  // also compute information matrix E(xx')

  int t;
  int i, j;
  int n = order[0];
  int m = order[1];
  int p = n + m + 1;
  int d = *delay;
  double *xt, *zt;

  xt = (double *)R_alloc(p, sizeof(double));
  zt = (double *)R_alloc(p, sizeof(double));

  // initialise outputs
  for (i = 0; i < p * p; i++)
    xz[i] = 0;
  for (i = 0; i < p; i++)
    xy[i] = 0;
  for (i = 0; i < p * p; i++)
    xx[i] = 0;

  for (t = *warmup; t < *len; t++) {
    // form the vectors x_t and z_t
    if (!my_isok(Y[t]))
      goto badsriv;
    for (j = 1; j <= n; j++) {
      i = j - 1;
      if (!my_isok(X[t - j]) || !my_isok(Y[t - j]))
        goto badsriv;
      xt[i] = X[t - j];
      zt[i] = Y[t - j];
    }
    for (j = 0; j <= m; j++) {
      i = n + j;
      if (!my_isok(U[t - (j + d)]))
        goto badsriv;
      xt[i] = zt[i] = U[t - (j + d)];
    }
    // add outer product to coefficient matrix xz
    for (i = 0; i < p; i++) {
      for (j = 0; j < p; j++) {
        xz[i + p * j] += xt[i] * zt[j];
      }
    }
    // add scalar product of x_t with Y_t to response vector xy
    for (i = 0; i < p; i++) {
      xy[i] += xt[i] * Y[t];
    }
    // add outer product to information matrix xx
    for (i = 0; i < p; i++) {
      for (j = 0; j < p; j++) {
        xx[i + p * j] += xt[i] * xt[j];
      }
    }
  badsriv:
    continue;
  }
}

// NOT USED
void sriv_infomat(double *U, double *X, int *len, int *warmup, int *order,
                  int *delay, double *xx) {
  // x: instrumental variable { Xf[k-1] ... Xf[k-n], Uf[k] ... Uf[k-m] }
  // form the information matrix E(xx')

  int t;
  int i, j;
  int n = order[0];
  int m = order[1];
  int p = n + m + 1;
  double *xt;

  xt = (double *)R_alloc(p, sizeof(double));

  // initialise outputs
  for (i = 0; i < p * p; i++)
    xx[i] = 0;

  for (t = *warmup; t < *len; t++) {
    // form the vector x_t
    for (j = 1; j <= n; j++) {
      i = j - 1;
      xt[i] = X[t - j];
    }
    for (j = 0; j <= m; j++) {
      i = n + j;
      xt[i] = U[t - (j + *delay)];
    }
    // add outer product to information matrix XX
    for (i = 0; i < p; i++) {
      for (j = 0; j < p; j++) {
        xx[i + p * j] += xt[i] * xt[j];
      }
    }
  }
}

/*
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

R_CMethodDef cMethods[] = {
        {"ar1_tv", &ar1_tv, 5, {REALSXP, REALSXP, INTSXP, REALSXP, REALSXP}},
        {NULL, NULL, 0}
};

void
R_init_hydromad(DllInfo *info)
{
        R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
*/
