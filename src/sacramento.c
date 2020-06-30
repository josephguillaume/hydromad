#include "sacramento.h"
#include <R.h>

/* Trivial wrapper function */
/* Based on code from University of Arizona MOSCEM project */
/* Simplified by Felix Andrews <felix@nfrac.org> 2010-02-01 */

void sma_sac(double *P, double *E, int *n, double *xpar, double *etmult,
             double *dt, double *U, double *uztwc_0, double *uzfwc_0,
             double *lztwc_0, double *lzfsc_0, double *lzfpc_0, double *adimc_0,
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
  }
}
