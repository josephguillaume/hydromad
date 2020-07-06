// Copyright (c) 2009 Felix Andrews <felix@nfrac.org>

#include <R_ext/Arith.h>

#include <R.h>

#ifndef min
#define min(a, b) ((a < b) ? (a) : (b))
#define max(a, b) ((a < b) ? (b) : (a))
#endif

#define my_isok(x) (!ISNA(x) & !ISNAN(x))

void swimp_core(double *overflow, int *n, double *alpha, double *beta,
                double *Emm, double *Pmm, double *Ksat, double *e, double *g,
                double *Hmax, double *Amax, double *porosity, double *drainage,
                double *undrainedV, double *V, double *A, double *H, double *M,
                double *Iw) {
  int t;
  double Hg, ETg, Mf;
  double pi = acos(-1);
  for (t = 1; t < *n; t++) {
    // work out water level in ground
    // relative to base of wetland, H==0
    Hg = *Hmax - (1 / *porosity) * M[t - 1];
    // infiltration rate (drainage) from wetland
    // depends on pressure; reference = 100mm
    Iw[t] = *Ksat * (H[t - 1] - Hg) / 100;
    // if discharging into wetland, adjust mass for porosity
    if (H[t - 1] < Hg)
      Iw[t] = Iw[t] * *porosity;
    // evapo-transpiration (equation from CMD model)
    Mf = M[t - 1] - Pmm[t];
    ETg = *e * Emm[t] * min(1, exp(2 * (1 - Mf / *g)));
    ETg = max(ETg, 0);
    // mass balance of water level in ground (CMD)
    M[t] = M[t - 1] + ETg - Pmm[t] - Iw[t] * (A[t - 1] / *Amax);
    M[t] = max(M[t], 0);
    // mass balance of volume in wetland
    V[t] = V[t - 1] + overflow[t] + (Pmm[t] - Emm[t]) * A[t - 1] -
           Iw[t] * A[t - 1];
    V[t] = max(V[t], 0);
    // drainage of volume above drainage level
    V[t] = V[t] - max(V[t] - *undrainedV, 0) * *drainage;
    // convert volume to water level, then area
    // NOTE: these formulae are designed to avoid numerical explosions!
    // e.g. ^(1/beta) when beta = 0.01
    A[t] = pow((2 / *beta + 1) * V[t] / (*alpha * pow(pi, -*beta / 2)),
               (2 / (2 + *beta)));
    H[t] = pow((2 / *beta + 1) * V[t] / pi, (*beta / (2 + *beta))) *
           pow(*alpha, (2 / (2 + *beta)));
  }
}
