/*
HBV model code for hydromad
R and Rcpp code by Alexander Buzacott (abuz5257@uni.sydney.edu.au)

Implementation based on HBV light as described in:
Seibert, J. and Vis, M. (2012). Teaching hydrological modeling with a user-
friendly catchment runoff-model software package. Hydrology and Earth System
Sciences, 16, 3315–3325, 2012.

HBV references:
Bergström, S.: The HBV Model: Its Structure and Applications,Swedish
Meteorological and Hydrological Institute (SMHI), Hydrology, Norrköping, 35
pp., 1992.

Bergström, S.: The HBV model (Chapter 13), in: Computer Models of Watershed
Hydrology, edited by: Singh, V. P., Water Resources Publications, Highlands
Ranch, Colorado, USA, 443–476, 1995.
*/

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame hbv_sim(NumericVector P, NumericVector E, NumericVector Tavg,
                  double tt, double cfmax, double sfcf, double cwh, double cfr,
                  double fc, double lp, double beta) {
  // Length of time series
  int nDays = P.size();

  // Storage output vector
  NumericVector U(nDays);

  // Set up vectors
  NumericVector sm = rep(0.0, nDays);
  NumericVector sp = rep(0.0, nDays);
  NumericVector wc = rep(0.0, nDays);
  NumericVector ETa = rep(0.0, nDays);
  NumericVector ep = rep(0.0, nDays);

  // Snow routine variables
  double infil;
  double refr;
  double maxwc;
  double melt;
  double snow;

  // Soil variables
  double recharge;
  double AET;
  double runoff;
  double smw;

  // Run model, starting at day 2
  for (int t = 1; t < nDays; t++) {
    // -----------------------------------------------------------------------
    // Snow routine
    // -----------------------------------------------------------------------
    infil = 0.0;

    // Determine if snow or rain falls
    if (Tavg[t] < tt) {
      // Refreezing of liquid in snow pack
      refr = std::min(cfr * cfmax * (Tavg[t] - tt), wc[t - 1]);
      wc[t] = wc[t - 1] - refr;
      // Use snowfall correction factor
      snow = P[t] * sfcf;
      // Add snow and refreezing water to snow pack
      sp[t] = sp[t - 1] + snow + refr;
      infil = P[t] - snow;
    } else {
      // Calculate and remove snowmelt
      melt = std::min(cfmax * (Tavg[t] - tt), sp[t - 1]);
      sp[t] = sp[t - 1] - melt;
      // Add melting water to liquid in snow pack
      wc[t] = wc[t - 1] + melt;
      // Calculate maximum liquid water holding capacity of snow pack
      maxwc = std::max(sp[t] * cwh, 0.0);
      if (wc[t] > maxwc) {
        // Add liquid excess water to effective P
        infil = P[t] + wc[t] - maxwc;
        wc[t] = maxwc;
      } else {
        // Retain liquid water in snow pack and just add P
        infil = P[t];
      }
    }
    // -----------------------------------------------------------------------
    // Soil routine
    // -----------------------------------------------------------------------
    recharge = 0.0;
    AET = 0.0;
    runoff = 0.0;

    // Calculate current soil wetness
    smw = std::max(std::min(pow(sm[t - 1] / fc, beta), 1.0), 0.0);

    // Calculate recharge and take away from infil
    recharge = infil * smw;
    infil -= recharge;

    // Add remaining infiltration to soil
    sm[t] = sm[t - 1] + infil;

    // Check if sm exceeds fc
    if (sm[t] >= fc) {
      runoff = sm[t] - fc;
      sm[t] = fc;
    }

    // Calculate actual ET
    AET = E[t] * std::min(sm[t] / (fc * lp), 1.0);
    if (AET < 0.0)
      AET = 0.0;
    if (sm[t] > AET) {
      ETa[t] = AET;
      sm[t] = sm[t] - AET;
    } else {
      ETa[t] = sm[t];
      sm[t] = 0.0;
    }
    ep[t] = runoff + recharge;
  }

  // Return
  return DataFrame::create(Named("U") = ep, Named("ETa") = ETa,
                           Named("sm") = sm, Named("sp") = sp,
                           Named("wc") = wc);
}

// [[Rcpp::export]]
DataFrame hbvrouting_sim(NumericVector U, double k0, double k1, double k2,
                         double uzl, double perc, NumericVector wi,
                         int n_maxbas) {
  // Length of timeseries
  int nDays = U.size();

  // Set up vectors
  NumericVector uz = rep(0.0, nDays);
  NumericVector lz = rep(0.0, nDays);
  NumericVector Qsim = rep(0.0, nDays);
  NumericVector X = rep(0.0, nDays);

  // Routing variables
  double actPERC;
  double Q0;
  double Q1;
  double Q2;

  for (int t = 1; t < nDays; t++) {
    // -----------------------------------------------------------------------
    // Discharge
    // -----------------------------------------------------------------------
    Q0 = 0.0;
    Q1 = 0.0;
    Q2 = 0.0;

    uz[t] = uz[t - 1] + U[t];

    // Percolation of water from upper to lower storage
    actPERC = std::min(uz[t], perc);
    uz[t] -= actPERC;
    lz[t] = lz[t - 1] + actPERC;

    // Quick runoff
    Q0 = k0 * std::max(uz[t] - uzl, 0.0);
    uz[t] -= Q0;

    Q1 = k1 * uz[t];
    uz[t] -= Q1;

    Q2 = k2 * lz[t];
    lz[t] -= Q2;

    Qsim[t] = Q0 + Q1 + Q2; // Total Q

    if (t < (n_maxbas - 1)) {
      X[t] = Qsim[t];
    } else {
      X[t] = sum(wi * (Qsim[Range(t - n_maxbas + 1, t)]));
    }
  }

  return DataFrame::create(Named("X") = X, Named("uz") = uz, Named("lz") = lz);
}