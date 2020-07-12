// hydromad: Hydrological Modelling and Analysis of Data
//
// Copyright (c) Felix Andrews <felix@nfrac.org>

#include <R_ext/Arith.h>

#include <R.h>
#include <R_ext/Print.h>/* for diagnostics */
#include <Rdefines.h>   /* need this? */

#ifndef min
#define min(a, b) ((a < b) ? (a) : (b))
#define max(a, b) ((a < b) ? (b) : (a))
#endif

#define my_isok(x) (!ISNA(x) & !ISNAN(x))

void inverse_filter(double *Q, int *nQ, double *a, double *b, int *n, int *m,
                    int *useQm, double *P, double *U) {
  int t, i;
  double Ut;
  for (t = max(*n, *m); t < *nQ; t++) {
    Ut = Q[t];
    for (i = 1; i <= *n; i++) {
      Ut -= a[i - 1] * Q[t - i];
    }
    for (i = 1; i <= *m; i++) {
      Ut -= b[i] * U[t - i];
    }
    /*
            if (*n >= 1) Ut = Ut - a[0] * Q[t-1];
            if (*n >= 2) Ut = Ut - a[1] * Q[t-2];
            if (*n >= 3) Ut = Ut - a[2] * Q[t-3];
            if (*m >= 1) Ut = Ut - b[1] * U[t-1];
            if (*m >= 2) Ut = Ut - b[2] * U[t-2];
            if (*m >= 3) Ut = Ut - b[3] * U[t-3];
    */
    Ut = Ut / b[0];
    U[t] = max(0, min(Ut, P[t]));
    if (*useQm > 0) {
      Q[t] = 0;
      for (i = 1; i <= *n; i++) {
        Q[t] += a[i - 1] * Q[t - i];
      }
      for (i = 0; i <= *m; i++) {
        Q[t] += b[i] * U[t - i];
      }
    }
  }
}

void inverse_filter_lambda(double *Q, int *nQ, double *Qs_0, double *Qq_0,
                           double *tau_s, double *tau_q, double *v_s_0,
                           double *lambda, double *P, double *U) {
  int t;
  double alpha_s = exp(-1 / *tau_s);
  double alpha_q = exp(-1 / *tau_q);
  double a = *v_s_0 * (alpha_q - alpha_s);
  double b = (1 - alpha_q);
  double c;
  double v_s_t, beta_s_t, beta_q_t;
  double Ubound0, Ubound1, Utest;
  double cbound0, cbound1, ctest;
  int positive_fn;
  double Qst1 = *Qs_0;
  double Qqt1 = *Qq_0;
  for (t = 0; t < *nQ; t++) {
    // trivial case: P is zero, therefore U is zero
    //	Rprintf("t = %d\n", t);
    if (P[t] == 0) {
      U[t] = 0;
      Qst1 = alpha_s * Qst1;
      Qqt1 = alpha_q * Qqt1;
      continue;
    }
    c = Q[t] - alpha_s * Qst1 - alpha_q * Qqt1;
    Ubound0 = c / (a + b);
    Ubound1 = (c - a) / b;
    Ubound0 = max(0, min(Ubound0, P[t]));
    Ubound1 = max(0, min(Ubound1, P[t]));
    // simple case: U is equal to P (both bounds exceed P)
    /*
            if ((Ubound0 >= P[t]) && (Ubound1 >= P[t])) {
                U[t] = P[t];
                if (U[t] > 0)
                    v_s_t = *v_s_0 * pow(U[t], *lambda);
                else v_s_t = 0;
                v_s_t = max(-1, min(0, v_s_t));
                beta_s_t = v_s_t * (1 - alpha_s);
                beta_q_t = (1 - v_s_t) * (1 - alpha_q);
                Qst1 = Qst1 + beta_s_t * U[t];
                Qqt1 = Qqt1 + beta_q_t * U[t];
                continue;
            }
    */

    // TODO: bound by U = -0.5
    // U[t] <- polyroot(-c, a, b)

    // hard case: need to search for U within possible bounds
    cbound0 = a * pow(Ubound0, 1.0 + *lambda) + b * Ubound0;
    cbound1 = a * pow(Ubound1, 1.0 + *lambda) + b * Ubound1;
    positive_fn = (cbound0 <= cbound1);
    Rprintf("U range: %.3f to %.3f\n", Ubound0, Ubound1);
    Rprintf("c range: %.3f to %.3f\n", cbound0, cbound1);
    // check if c is outside possible range, just take closest value
    /*	if ((c <= cbound0) && (c <= cbound1)) { // take min of cbounds
                if (positive_fn) Ubound1 = Ubound0;
                else Ubound0 = Ubound1;
            }
            if ((c >= cbound0) && (c >= cbound1)) { // take max of cbounds
                if (positive_fn) Ubound0 = Ubound1;
                else Ubound1 = Ubound0;
            }
    */
    // otherwise, do iterative bisection of U range to match c
    while (fabs(Ubound0 - Ubound1) > 0.0001) {
      // choose new U value between existing bounds
      // TODO: use linear approximation rather than bisection
      Utest = (Ubound0 + Ubound1) * 0.5;
      ctest = a * pow(Utest, 1.0 + *lambda) + b * Utest;
      Rprintf("U optim loop: Ut = %.3f, ct = %.3f, c = %.3f\n", Utest, ctest,
              c);
      if ((positive_fn && (ctest > c)) || (!positive_fn && (ctest < c))) {
        cbound1 = ctest;
        Ubound1 = Utest;
      } else {
        cbound0 = ctest;
        Ubound0 = Utest;
      }
    }
    U[t] = Ubound0;
    v_s_t = 0.0;
    if (U[t] > 0)
      v_s_t = *v_s_0 * pow(U[t], *lambda);
    v_s_t = max(0.0, min(1.0, v_s_t));
    beta_s_t = v_s_t * (1.0 - alpha_s);
    beta_q_t = (1.0 - v_s_t) * (1.0 - alpha_q);
    Qst1 = alpha_s * Qst1 + beta_s_t * U[t];
    Qqt1 = alpha_q * Qqt1 + beta_q_t * U[t];
  }
}

double ctestfun(double a, double b, double U, double lambda) {
  return (a * pow(U, 1 + lambda) + b * U);
}
