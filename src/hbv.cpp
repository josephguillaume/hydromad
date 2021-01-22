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
                  double tt, double cfmax, double sfcf, double cfr, double cwh,
                  double fc, double lp, double beta, bool initialise_sm) {
  // Length of time series
  int nDays = P.size();

  // Storage output vector
  NumericVector U(nDays);

  // Set up vectors and variables
  NumericVector sm = rep(0.0, nDays);
  NumericVector sp = rep(0.0, nDays);
  NumericVector infil = rep(0.0, nDays);
  NumericVector AET = rep(0.0, nDays);
  NumericVector recharge = rep(0.0, nDays);

  // Snow routine variables
  double refr = 0.0;
  double wc_ = 0.0;
  double sp_ = 0.0;
  double melt;

  double infil_;
  double sp_tm1;
  double maxwc;

  // Soil variables
  double sm_ = 0.0;
  if (initialise_sm)
    sm_ = fc * lp;
  double sm_tm1;
  double infil_r;
  double infil_s;
  double rm;
  double sm_et;

  // Run model
  for (int t = 0; t < nDays; t++) {
    // -----------------------------------------------------------------------
    // Snow routine
    // -----------------------------------------------------------------------
    infil_ = 0.0;
    sp_tm1 = sp_;
    if (P[t] > 0) {
      // Determine if snow or rain falls
      if (Tavg[t] > tt) {
        wc_ += P[t];
      } else {
        // Snow and apply snowfall correction factor
        sp_ += P[t] * sfcf;
      }
    }

    if (Tavg[t] > tt) {
      // Melt snow
      melt = cfmax * (Tavg[t] - tt);
      // If melt is greater than snow depth
      if (melt > sp_) {
        // All water is added to infiltration
        infil_ = sp_ + wc_;
        wc_ = 0;
        sp_ = 0;
      } else {
        // Remove melt from snow pack
        sp_ -= melt;
        wc_ += melt;
        // Calculate maximum liquid water holding capacity of snow pack
        maxwc = sp_ * cwh;
        if (wc_ > maxwc) {
          infil_ = wc_ - maxwc;
          wc_ = maxwc;
        }
      }
    } else {
      // Refreeze water in liquid snow store
      refr = std::min(cfr * cfmax * (tt - Tavg[t]), wc_);
      sp_ += refr;
      wc_ -= refr;
    }
    sp[t] = sp_ + wc_;

    // -----------------------------------------------------------------------
    // Soil routine
    // -----------------------------------------------------------------------
    // Divide portion of infiltration that goes to soil/gw
    sm_tm1 = sm_;
    if (infil_ > 0) {
      if (infil_ < 1) {
        infil_s = infil_;
      } else {
        infil_r = std::round(infil_);
        infil_s = infil_ - infil_r;
        int i = 1;
        while (i <= infil_r) {
          rm = std::pow(sm_ / fc, beta);
          if (rm > 1)
            rm = 1;
          sm_ += 1 - rm;
          recharge[t] += rm;
          i++;
        }
      }
      rm = std::pow(sm_ / fc, beta);
      if (rm > 1)
        rm = 1;
      sm_ += (1 - rm) * infil_s;
      recharge[t] += rm * infil_s;
    }
    // Only AET if there is snow cover the previous timestep
    if (sp_tm1 == 0) {
      sm_et = (sm_ + sm_tm1) / 2;
      // Calculate actual ET
      AET[t] = E[t] * std::min(sm_et / (fc * lp), 1.0);
      if (AET[t] < 0)
        AET[t] = 0;
      // Remove AET from soil if there is water
      if (sm_ > AET[t]) {
        sm_ -= AET[t];
      } else {
        AET[t] = sm_;
        sm_ = 0;
      }
    }
    sm[t] = sm_;
  }

  // Return
  return DataFrame::create(Named("U") = recharge, Named("sp") = sp,
                           Named("sm") = sm, Named("AET") = AET);
}

// [[Rcpp::export]]
DataFrame hbvrouting_sim(NumericVector U, double perc, double uzl, double k0,
                         double k1, double k2, NumericVector wi, int n_maxbas,
                         double initial_lz) {
  // Length of timeseries
  int nDays = U.size();

  // Set up vectors
  NumericVector suz = rep(0.0, nDays);
  NumericVector slz = rep(0.0, nDays);
  NumericVector Q0 = rep(0.0, nDays);
  NumericVector Q1 = rep(0.0, nDays);
  NumericVector Q2 = rep(0.0, nDays);

  // Routing variables
  double act_perc = 0;
  double suz_ = 0;
  double slz_ = initial_lz;

  for (int t = 0; t < nDays; t++) {
    // -----------------------------------------------------------------------
    // Discharge
    // -----------------------------------------------------------------------
    // Add runoff and recharge to upper zone of storage
    suz_ += U[t];
    // Percolation of of water from upper to lower zone
    act_perc = std::min(suz_, perc);
    suz_ -= act_perc;
    slz_ += act_perc;

    // Calculate runoff from storage
    Q0[t] = k0 * std::max(suz_ - uzl, 0.0);
    Q1[t] = k1 * suz_;
    suz_ -= (Q1[t] + Q0[t]);

    Q2[t] = k2 * slz_;
    slz_ -= Q2[t];

    suz[t] = suz_;
    slz[t] = slz_;
  }

  return DataFrame::create(Named("suz") = suz, Named("slz") = slz,
                           Named("Q0") = Q0, Named("Q1") = Q1,
                           Named("Q2") = Q2);
}

// [[Rcpp::export]]
NumericVector hbv_pet(DateVector dates, NumericVector Tavg, NumericVector PET,
                      NumericVector Tmean, double cet) {
  // Setup variables
  int nDays = Tavg.size();
  int doy;
  int year;
  double pet_;
  Date d;
  NumericVector E(nDays);

  for (int i = 0; i < nDays; i++) {
    d = dates[i];
    doy = d.getYearday() - 1;
    year = d.getYear();
    if ((year % 4 == 0) & ((year % 100 != 0) | (year % 400 == 0))) {
      if (doy > 58)
        doy -= 1;
    }
    pet_ = 1.0 + cet * (Tavg[i] - Tmean[doy]);
    if (pet_ > 2.0)
      pet_ = 2.0;
    if (pet_ < 0.0)
      pet_ = 0.0;
    E[i] = pet_ * PET[doy];
  }
  return E;
}