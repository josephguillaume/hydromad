#include <Rcpp.h>
using namespace Rcpp;

// based on Chiew et al. 2009

// [[Rcpp::export]]
DataFrame simhyd_sim(NumericVector P, NumericVector E, double INSC,
                     double COEFF, double SQ, double SMSC, double SUB,
                     double CRAK, double K, double GWt0, double SMSt0) {
  // size
  int n = P.size();
  // storage output vector
  NumericVector U(n);
  // component flux vectors
  NumericVector IMAX(n), INT(n), INR(n), RMO(n), IRUN(n);
  NumericVector ET(n), SRUN(n), REC(n), SMF(n), POT(n), BAS(n);
  // Stores
  NumericVector SMS(n), GW(n);
  // at previous timestep
  double SMSt1, GWt1;

  // initialise vectors
  SMSt1 = SMSC * SMSt0;
  GWt1 = GWt0;
  // Rcout << SMSt1;

  // run a loop
  for (int t = 0; t < n; ++t) {
    // Interception storage
    IMAX[t] = std::min(INSC, E[t]);
    // Interception
    INT[t] = std::min(IMAX[t], P[t]);
    // interception runoff (INR)
    INR[t] = P[t] - INT[t];
    // Infiltration capacity
    RMO[t] = std::min(COEFF * std::exp(-SQ * SMSt1 / SMSC), INR[t]);
    // Rcout << RMO[t];
    // Direct runoff
    IRUN[t] = INR[t] - RMO[t];
    // Saturation Excess runoff (SRUN)
    SRUN[t] = SUB * SMSt1 / SMSC * RMO[t];
    // Recharge
    REC[t] = CRAK * SMSt1 / SMSC * (RMO[t] - SRUN[t]);
    // Infiltration into soil store
    SMF[t] = RMO[t] - SRUN[t] - REC[t];
    // potential ET
    POT[t] = E[t] - INT[t];
    // calculate soil ET
    ET[t] = std::min(10 * SMSt1 / SMSC, POT[t]);
    // calulate SMS overflow (see figure 2 Chiew et al. 2009)
    // calculate soil moisture balance
    SMS[t] = SMSt1 + SMF[t] - ET[t];
    if (SMS[t] > SMSC) {
      SMS[t] = SMSC;
      REC[t] = REC[t] + SMS[t] - SMSC;
    }
    SMSt1 = SMS[t];
    // 		Rcout << SMS[t];
    // calculate baseflow
    BAS[t] = K * GWt1;
    // calculate GW storage
    GW[t] = GWt1 + REC[t] - BAS[t];
    GWt1 = GW[t];
    // Calculate total local runoff
    U[t] = IRUN[t] + SRUN[t] + BAS[t];
    // 		Rcout << U[t];
  }
  // combine with ET in a matrix
  return DataFrame::create(Named("U") = U, Named("ET") = ET);
}
