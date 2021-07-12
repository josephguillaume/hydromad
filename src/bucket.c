// hydromad: Hydrological Modelling and Analysis of Data
//
// Copyright (c) Felix Andrews <felix@nfrac.org>

#include <R_ext/Arith.h>

#include <R.h>

#ifndef min
#define min(a, b) ((a < b) ? (a) : (b))
#define max(a, b) ((a < b) ? (b) : (a))
#endif

void sma_bucket(double *P, double *E, int *n, double *Sb, double *fc,
                double *Sfc, double *a_ei, double *M, double *a_ss, double *S_0,
                double *U, double *S, double *ET) {
  int t;
  double Eintc, Etrans, Ebare, Use, Uss;
  double S_prev = *S_0;
  for (t = 0; t < *n; t++) {
    // evapo-transpiration
    Eintc = *a_ei * P[t];
    S[t] = min(*Sb, S_prev + P[t] - Eintc);
    Etrans = *M * min(1, S[t] / *Sfc) * E[t];
    Ebare = (1 - *M) * (S[t] / *Sb) * E[t];
    ET[t] = Eintc + min(S[t], Etrans + Ebare);
    // mass balance
    S[t] = S_prev + P[t] - ET[t];
    // drainage (saturation excess)
    Use = max(0, S[t] - *Sb);
    S[t] = S[t] - Use;
    // drainage (sub-surface)
    Uss = max(0, *a_ss * (S[t] - *Sfc));
    S[t] = S[t] - Uss;
    U[t] = Use + Uss;
    S_prev = S[t];
  }
}
