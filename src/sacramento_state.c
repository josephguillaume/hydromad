#include "sacramento_state.h"
#include <R.h>

/* Trivial wrapper function */
/* Based on code from University of Arizona MOSCEM project */
/* Simplified by Felix Andrews <felix@nfrac.org> 2010-02-01 */
/* Adapted to return state by Joseph Guillaume 2013-10-24 */

void sma_sac_state(double *P, double *E, int *n, double *xpar, double *etmult,
                   double *dt, double *U, double *uztwc, double *uzfwc,
                   double *lztwc, double *lzfsc, double *lzfpc, double *adimc,
                   double *sett, double *se1, double *se3, double *se4,
                   double *se5, double *roimp, double *sdro, double *ssur,
                   double *sif, double *bfp, double *bfs, double *bfcc,
                   double *uztwc_0, double *uzfwc_0, double *lztwc_0,
                   double *lzfsc_0, double *lzfpc_0, double *adimc_0,
                   int *min_ninc) {
  struct SMA sma;
  struct FSUM1 fsum1;

  /* ASSIGN MODEL SPECIFIC PARAMETERS TO PARAMETERS STRUCT */
  sma.uztwm = xpar[0];
  sma.uzfwm = xpar[1];
  sma.uzk = xpar[2];
  sma.pctim = xpar[3];
  sma.adimp = xpar[4];
  sma.riva = 0.0;
  sma.zperc = xpar[5];
  sma.rexp = xpar[6];
  sma.lztwm = xpar[7];
  sma.lzfsm = xpar[8];
  sma.lzfpm = xpar[9];
  sma.lzsk = xpar[10];
  sma.lzpk = xpar[11];
  sma.pfree = xpar[12];
  sma.rserv = 0.3;
  sma.side = 0.0;
  sma.pxmlt = 1.0;

  sma.uztwc = *uztwc_0 * sma.uztwm;
  sma.uzfwc = *uzfwc_0 * sma.uzfwm;
  sma.lztwc = *lztwc_0 * sma.lztwm;
  sma.lzfsc = *lzfsc_0 * sma.lzfsm;
  sma.lzfpc = *lzfpc_0 * sma.lzfpm;
  sma.adimc = *adimc_0 * (sma.uztwm + sma.lztwm);

  sma.min_ninc = *min_ninc;

  /* SET SOME INITIAL VALUES TO ZERO */
  fsum1.srot = fsum1.simpvt = fsum1.srodt = fsum1.srost = 0.;
  fsum1.sintft = fsum1.sgwfp = fsum1.sgwfs = fsum1.srecht = fsum1.sett = 0.;
  fsum1.se1 = fsum1.se3 = fsum1.se4 = fsum1.se5 = 0.;

  /* SET EVAPOTRANSPIRATION DISTRIBUTION */
  sma.epdist = *etmult;

  /* DT IS THE LENGTH OF EACH TIME INTERVAL IN DAYS */
  /* DT IS USED TO CALCULATE dinc IN fland1 */
  sma.dt = *dt; /*  */

  for (int t = 0; t < *n; t++) {

    /* ASSIGN AND ADJUST PET VALUE */
    sma.ep = E[t] * sma.epdist;

    /* ASSIGN AND ADJUST PRECIPITATION VALUE */
    sma.pxv = P[t] * sma.pxmlt;

    /* PERFORM SOIL MOISTURE ACCOUNTING OPERATIONS  */
    fland1(&sma, &fsum1);

    /* SET TOTAL CHANNEL INFLOW EQUAL TO THE EXCESS AT THE END  */
    /* OF EACH TIME PERIOD */
    U[t] = sma.tlci;

    /* SAVE STATE VARIABLES */
    uztwc[t] = sma.uztwc;
    uzfwc[t] = sma.uzfwc;
    lztwc[t] = sma.lztwc;
    lzfsc[t] = sma.lzfsc;
    lzfpc[t] = sma.lzfpc;
    adimc[t] = sma.adimc;

    /* SAVE EVAPORATION TOTALS */
    sett[t] = fsum1.sett;
    se1[t] = fsum1.se1;
    se3[t] = fsum1.se3;
    se4[t] = fsum1.se4;
    se5[t] = fsum1.se5;

    /* SAVE FLOWS */
    roimp[t] = sma.tlci_flows[0];
    sdro[t] = sma.tlci_flows[1];
    ssur[t] = sma.tlci_flows[2];
    sif[t] = sma.tlci_flows[3];
    bfp[t] = sma.tlci_flows[4];
    bfs[t] = sma.tlci_flows[5];
    bfcc[t] = sma.tlci_flows[6];
  }
}
