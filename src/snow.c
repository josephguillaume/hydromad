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

// converted by Felix from code by Jarkko Koskela @tkk.fi 2010-02-26

void sma_snow(double *Prain, double *Psnow, double *E, int *n, double *kd,
              double *kf, double *rcap, double *Tmelt, double *LSWE_0,
              double *ISWE_0, double *U, double *SWE) {
  int t;
  double melt, freeze, ISWE, LSWE;
  double LSWEprev = *LSWE_0;
  double ISWEprev = *ISWE_0;
  for (t = 0; t < *n; t++) {
    // Melt (degree day model)
    melt = min(max(*kd * (E[t] - *Tmelt), 0), ISWEprev);
    // Freezing (degree day model)
    freeze = min(max(*kf * (*Tmelt - E[t]), 0), LSWEprev);
    // Mass balance for the snowpack
    //
    // Ice in the snowpack
    ISWE = ISWEprev + Psnow[t] + freeze - melt;
    // Water in the snowpack
    LSWE = min(*rcap * ISWE, LSWEprev + Prain[t] + melt - freeze);
    // Rain/melt is snowmelt discharge when there is snow on the ground,
    // and rainfall in snow-free periods.
    U[t] = max(Prain[t] + melt - freeze - (*rcap * ISWE - LSWEprev), 0);
    SWE[t] = LSWE + ISWE;
    ISWEprev = ISWE;
    LSWEprev = LSWE;
  }
}
