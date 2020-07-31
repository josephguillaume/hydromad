#include "sacramento.h"
#include <R.h>
#include <R_ext/Print.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define EPS 1.e-5 /* Minimum parameter value */

void fland1(struct SMA *sma, struct FSUM1 *fsum1) {

  /* THIS SUBROUTINE EXECUTES THE SAC-SMA OPERATION FOR ONE
   * TIME PERIOD.
   ********************************************************
   * SUBROUTINE INITIALLY WRITTEN BY. . .
   * ERIC ANDERSON - HRL     APRIL 1979     VERSION 1
   * MODIFIED FEB 1988 FOR RESEARCH ON THE U OF A CAMPUS */
  /*********************************************************/

  /*  WRITTEN IN C. FROM FORTRAN CODE BY
   *  PATRICE O. YAPO 7/27/93  */
  /*  Last revised poy 2/20/94 */

  int i, ninc;
  double edmnd, e1, e2, e3, e4, e5, red, uzrat, ratlzt, saved, ratlz;
  double del, twx, roimp, ssur, sif, sperc, sdro, spbf, sbf;
  double dinc, pinc, duz, dlzp, dlzs, parea, adsur, ratio;
  double addro, bf, percm, zp, perc, defr, check, perct;
  // double uzdefr; //upper zone deficiency ratio, not used
  double percf, hpl, ratlp, ratls, fracp, percp, percs, excess;
  double sur, eused, tbf, bfcc, bfp, bfs, bfncc, tet;
  double rsum[7];

  /*    REPARAMETERIZATION OPTION */
  /* izflag=2; reparameterized percolation,
   *           otherwise, original percolation */
  /* izflag = 1;
     if (izflag == 2) aperc = sma->zperc;
  */

  /*   if (iter1>=5050) printf("enter fland1\n");*/

  zp = sma->zperc;

  /* COMPUTE EVAPOTRANSPIRATION LOSS FOR THE TIME INTERVAL.
   * EDMND IS THE ET-DEMAND FOR THE TIME INTERVAL */

  /* edmnd = sma->ep * sma->epdist; */ /* Compute et from upper zone */
  edmnd = sma->ep;
  e1 = edmnd * sma->uztwc / sma->uztwm;
  red = edmnd - e1; /* red is residual evap demand */
  sma->uztwc -= e1;
  if (fabs(sma->uztwc) < EPS)
    sma->uztwc = 0.;
  e2 = 0.;
  if (sma->uztwc < 0.) { /* e1 cannot exceed sma->uztwc */
    e1 += sma->uztwc;
    sma->uztwc = 0.;
    red = edmnd - e1;
    ;
    if (sma->uzfwc < red) {
      e2 = sma->uzfwc; /* e2 is evap from sma->uzfwc */
      sma->uzfwc = 0.;
      red -= e2;
      goto L225;
    } else {
      e2 = red;
      sma->uzfwc -= e2;
      red = 0.;
    }
  }

  /* UPPER ZONE FREE WATER RATIO EXCEEDS UPPER ZONE TENSION
   * WATER RATIO, THUS TRANSFER FREE WATER TO TENSION */

  if ((sma->uztwc / sma->uztwm) < (sma->uzfwc / sma->uzfwm)) {
    uzrat = (sma->uztwc + sma->uzfwc) / (sma->uztwm + sma->uzfwm);
    sma->uztwc = sma->uztwm * uzrat;
    sma->uzfwc = sma->uzfwm * uzrat;
  }
L225: /* This is a C label */
  e3 = red * sma->lztwc / (sma->uztwm + sma->lztwm);
  sma->lztwc -= e3;
  if (fabs(sma->lztwc) < EPS)
    sma->lztwc = 0.;
  if (sma->lztwc < 0.) {
    e3 += sma->lztwc; /* e3 cannot exceed sma->lztwc */
    sma->lztwc = 0.;
  }
  ratlzt = sma->lztwc / sma->lztwm;
  saved = sma->rserv * (sma->lzfpm + sma->lzfsm);
  ratlz = (sma->lztwc + sma->lzfpc + sma->lzfsc - saved) /
          (sma->lztwm + sma->lzfpm + sma->lzfsm - saved);

  if (ratlzt < ratlz) { /* RESUPPLY LOWER ZONE TENSION
                          WATER FROM LOWER ZONE FREE WATER
                          IF MORE WATER AVAILABLE THERE. */
    del = (ratlz - ratlzt) * sma->lztwm;
    /*    TRANSFER FROM sma->lzfsc to sma->lztwc */
    sma->lztwc = sma->lztwc + del;

    sma->lzfsc = sma->lzfsc - del;
    if (sma->lzfsc < 0.) { /* IF TRANSFER EXCEEDS LZFSC THEN
                         REMAINDER COMES FROM LZFPC */
      sma->lzfpc += sma->lzfsc;
      sma->lzfsc = 0.;
    }
  }

  /*   COMPUTE ET FROM ADIMP AREA.-E5 */

  e5 = e1 +
       (red + e2) * (sma->adimc - e1 - sma->uztwc) / (sma->uztwm + sma->lztwm);

  /*    ADJUST ADIMC,ADDITIONAL IMPERVIOUS AREA STORAGE, FOR
   *    EVAPORATION */

  sma->adimc -= e5;
  if (fabs(sma->adimc) < EPS)
    sma->adimc = 0.;
  if (sma->adimc < 0.) { /* E5 CAN NOT EXCEED ADIMC. */
    e5 += sma->adimc;
    sma->adimc = 0.;
  }
  e5 *= sma->adimp;

  /*    E5 IS ET FROM THE AREA ADIMP.
   *    COMPUTE PERCOLATION AND RUNOFF AMOUNTS. */

  twx = sma->pxv + sma->uztwc - sma->uztwm;
  if (twx < 0.) { /* ALL MOISTURE HELD IN UZTW--NO EXCESS.  */
    sma->uztwc += sma->pxv;
    twx = 0.;
  } else { /* MOISTURE AVAILABLE IN EXCESS OF UZTW STORAGE */
    sma->uztwc = sma->uztwm;
  }
  sma->adimc += (sma->pxv - twx);

  /*    COMPUTE IMPERVIOUS AREA RUNOFF. */
  roimp = sma->pxv * sma->pctim;

  /*    ROIMP IS RUNOFF FROM THE MINIMUM IMPERVIOUS AREA.  */
  fsum1->simpvt += roimp;

  /*    INITIALIZE TIME INTERVAL SUMS. */
  sbf = ssur = sif = sperc = sdro = spbf = 0.0;

  /*    DETERMINE COMPUTATIONAL TIME INCREMENTS FOR THE BASIC
   *    TIME INTERVAL */
  /*    NINC = NUMBER OF TIME INCREMENTS THAT THE TIME INTERVAL
   *    IS DIVIDED INTO FOR FURTHER SOIL-MOISTURE ACCOUNTING.
   *    NO ONE INCREMENT WILL EXCEED 5.0 MILLIMETERS OF UZFWC+PAV
   *    DINC = LENGTH OF EACH INCREMENT IN DAYS.
   *    PINC = AMOUNT OF AVAILABLE MOISTURE FOR EACH INCREMENT.
   *     ninc = (int) (1.0 + 0.2 * sma->uzfwc + twx);      */

  ninc = (int)(1.0 + 0.2 * (sma->uzfwc + twx));
  if (ninc < sma->min_ninc)
    ninc = sma->min_ninc; // minimum number of loops
  dinc = 1.0 / (double)ninc * sma->dt;
  pinc = twx / (double)ninc;

  /*    COMPUTE FREE WATER DEPLETION FRACTIONS FOR
   *    THE TIME INCREMENT BEING USED-BASIC DEPLETIONS
   *    ARE FOR ONE DAY */

  if (sma->uzk > 1.)
    Rprintf("uzk = %f\n", sma->uzk);
  if (sma->lzpk > 1.)
    Rprintf("lzpk = %f\n", sma->lzpk);
  if (sma->lzsk > 1.)
    Rprintf("lzsk = %f\n", sma->lzsk);
  duz = 1.0 - pow((1.0 - sma->uzk), dinc);
  dlzp = 1.0 - pow((1.0 - sma->lzpk), dinc);
  dlzs = 1.0 - pow((1.0 - sma->lzsk), dinc);
  parea = 1.0 - sma->adimp - sma->pctim;

  /*    START INCREMENTAL FOR LOOP FOR THE TIME INTERVAL.  */
  /*      if (iter1>=5050) printf("ninc=%d\n",ninc);*/

  for (i = 0; i < ninc; i++) {
    adsur = 0.;

    /*        COMPUTE DIRECT RUNOFF (FROM ADIMP AREA).
     *        ADDRO IS THE AMOUNT OF DIRECT RUNOFF FROM
     *        THE AREA ADIMP-SDRO IS THE SIX HOUR SUMMATION */
    ratio = (sma->adimc - sma->uztwc) / sma->lztwm;
    addro = pinc * (ratio * ratio);
    sdro += (addro * sma->adimp);

    /*        COMPUTE BASEFLOW AND KEEP TRACK OF TIME INTERVAL SUM.  */
    bf = sma->lzfpc * dlzp;
    sma->lzfpc -= bf;
    if (sma->lzfpc <= 1.e-4) {
      bf += sma->lzfpc;
      sma->lzfpc = 0.0;
    }
    sbf += bf;
    spbf += bf;
    bf = sma->lzfsc * dlzs;
    sma->lzfsc -= bf;
    if (sma->lzfsc <= 1.e-4) {
      bf += sma->lzfsc;
      sma->lzfsc = 0.0;
    }
    sbf += bf;

    /*        COMPUTE PERCOLATION-IF NO WATER AVAILABLE THEN SKIP */
    if ((pinc + sma->uzfwc) <= 1.e-2) {
      sma->uzfwc += pinc;
      goto L249;
    }
    percm = sma->lzfpm * dlzp + sma->lzfsm * dlzs;
    /*
              if (izflag == 2) {
                 zp = (sma->uzfwm - percm) / (percm * pow(aperc,sma->rexp));
              }
    */
    if (zp < 0.0)
      zp = 0.0;
    perc = percm * sma->uzfwc / sma->uzfwm;
    /*        DEFR IS THE LOWER ZONE MOISTURE DEFICIENCY RATIO */
    defr = 1.0 - (sma->lztwc + sma->lzfpc + sma->lzfsc) /
                     (sma->lztwm + sma->lzfpm + sma->lzfsm);
    if (defr < 0.) {
      Rprintf("defr = %f\n", defr);
      Rprintf("%f  %f  %f\n", sma->lztwc, sma->lzfpc, sma->lzfsc);
      Rprintf("%f  %f  %f\n", sma->lztwm, sma->lzfpm, sma->lzfsm);
      error("error: lower zone moisture deficiency ratio is < 0.");
    }
    // uzdefr = 1.0 - (sma->uztwc + sma->uzfwc) / (sma->uztwm + sma->uzfwm);
    perc *= (1.0 + zp * pow(defr, sma->rexp));

    /* NOTE...PERCOLATION OCCURS FROM UZFWC BEFORE PAV IS ADDED */

    if (perc >= sma->uzfwc) { /* PERCOLATION RATE EXCEEDS UZFWC */
      perc = sma->uzfwc;
    }
    sma->uzfwc -= perc; /* PERCOLATION RATE IS LESS THAN UZFWC.  */

    check = sma->lztwc + sma->lzfpc + sma->lzfsc + perc - sma->lztwm -
            sma->lzfpm - sma->lzfsm;
    if (check > 0.) { /* CHECK TO SEE IF PERCOLATION
                       * EXCEEDS LOWER ZONE DEFICIENCY.  */
      perc -= check;
      sma->uzfwc += check;
    }
    /*        SPERC IS THE TIME INTERVAL SUMMATION OF PERC */
    sperc += perc;

    /*        COMPUTE INTERFLOW AND KEEP TRACK OF TIME INTERVAL SUM.
     * NOTE...PINC HAS NOT YET BEEN ADDED */
    del = sma->uzfwc * duz;
    sif += del;
    sma->uzfwc -= del;

    /*        DESCRIBE PERCOLATED WATER INTO THE LOWER ZONES
     *        TENSION WATER MUST BE FILLED FIRST EXCEPT FOR THE
     *        PFREE AREA.  PERCT IS PERCOLATION TO TENSION WATER
     *        AND PERCF IS PERCOLATION GOING TO FREE WATER.  */
    perct = perc * (1.0 - sma->pfree);
    if ((perct + sma->lztwc) <= sma->lztwm) {
      sma->lztwc += perct;
      percf = 0.0;
    } else {
      percf = perct + sma->lztwc - sma->lztwm;
      sma->lztwc = sma->lztwm;
    }

    /*        DISTRIBUTE PERCOLATION IN EXCESS OF TENSION
     *        REQUIREMENTS AMONG THE FREE WATER STORAGES */
    percf += (perc * sma->pfree);
    if (percf != 0.) {

      /*           HPL IS THE RELATIVE SIZE OF THE PRIMARY STORAGE
       *           AS COMPARED WITH TOTAL LOWER ZONE FREE WATER STORAGE.  */
      hpl = sma->lzfpm / (sma->lzfpm + sma->lzfsm);

      /*           RATLP AND RATLS ARE CONTENT TO CAPACITY RATIOS, OR
       *           IN OTHER WORDS, THE RELATIVE FULLNESS OF EACH STORAGE */
      ratlp = sma->lzfpc / sma->lzfpm;
      ratls = sma->lzfsc / sma->lzfsm;
      /*           FRACP IS THE FRACTION GOING TO PRIMARY. */
      fracp = hpl * 2.0 * (1.0 - ratlp) / (1.0 - ratlp + 1.0 - ratls);
      if (fracp > 1.0)
        fracp = 1.0;

      /*           PERCP AND PERCS ARE THE AMOUNT OF THE EXCESS
       *           PERCOLATION GOING TO PRIMARY AND SUPPLEMENTAL
       *           STORAGES, RESPECTIVELY. */
      percp = percf * fracp;
      percs = percf - percp;
      sma->lzfsc += percs;
      if (sma->lzfsc > sma->lzfsm) {
        percs += (-sma->lzfsc + sma->lzfsm);
        sma->lzfsc = sma->lzfsm;
      }
      sma->lzfpc += (percf - percs);

      /*           CHECK TO MAKE SURE LZFPC DOES NOT EXCEED LZFPM */
      if (sma->lzfpc > sma->lzfpm) {
        excess = sma->lzfpc - sma->lzfpm;
        sma->lztwc += excess;
        sma->lzfpc = sma->lzfpm;
      }
    }

    /*        DISTRIBUTE PINC BETWEEN UZFWC AND SURFACE RUNOFF. */

    if (pinc != 0.) {
      if ((pinc + sma->uzfwc) <= sma->uzfwm) { /* CHECK IF PINC
                                                * EXCEEDS UZFWM */
        sma->uzfwc += pinc;                    /* NO SUFACE RUNOFF */
      } else {
        sur = pinc + sma->uzfwc - sma->uzfwm;
        sma->uzfwc = sma->uzfwm;

        /*              ADSUR IS THE AMOUNT OF SURFACE RUNOFF WHICH COMES
         *              FROM THAT PORTION OF ADIMP WHICH IS NOT
         *              CURRENTLY GENERATING DIRECT RUNOFF.  ADDRO/PINC
         *              IS THE FRACTION OF ADIMP CURRENTLY GENERATING
         *              DIRECT RUNOFF. */
        ssur = ssur + (sur * parea);
        adsur = sur * (1.0 - addro / pinc);
        ssur += (adsur * sma->adimp);
      }
    }
  L249: /* C Label */
    sma->adimc += (pinc - addro - adsur);
    if (sma->adimc > (sma->uztwm + sma->lztwm)) {
      addro += (sma->adimc - (sma->uztwm + sma->lztwm));
      sma->adimc = sma->uztwm + sma->lztwm;
    }
  } /* END OF INCREMENTAL FOR LOOP. */

  /*   COMPUTE SUMS AND ADJUST RUNOFF AMOUNTS BY THE AREA OVER
   *   WHICH THEY ARE GENERATED. */

  /*    EUSED IS THE ET FROM PAREA WHICH IS 1.0-ADIMP-PCTIM */
  eused = e1 + e2 + e3;
  sif *= parea;

  /*    SEPARATE CHANNEL COMPONENT OF BASEFLOW FROM THE
   *    NON-CHANNEL COMPONENT */

  tbf = sbf * parea; /* TBF IS TOTAL BASEFLOW */
                     /*    BFCC IS BASEFLOW, CHANNEL COMPONENT */
  bfcc = tbf * (1.0 / (1.0 + sma->side));
  bfp = (spbf * parea) / (1.0 + sma->side);
  bfs = bfcc - bfp;
  if (bfs < 0.)
    bfs = 0;
  bfncc = tbf - bfcc; /* BFNCC IS BASEFLOW, NON-CHANNEL COMPONENT */

  /*    ADD TO MONTHLY SUMS. */
  fsum1->sintft += sif;
  fsum1->sgwfp += bfp;
  fsum1->sgwfs += bfs;
  fsum1->srecht += bfncc;
  fsum1->srost += ssur;
  fsum1->srodt += sdro;

  /*    STORE EACH OF THE FIVE FLOWS IN VARIABLE TLCI_FLOWS*/ /* ADDED BY DPD */
  sma->tlci_flows[0] = roimp;                                 /* ADDED BY DPD */
  sma->tlci_flows[1] = sdro;                                  /* ADDED BY DPD */
  sma->tlci_flows[2] = ssur;                                  /* ADDED BY DPD */
  sma->tlci_flows[3] = sif;                                   /* ADDED BY DPD */
  sma->tlci_flows[4] = bfp;                                   /* ADDED BY DPD */
  sma->tlci_flows[5] = bfs;                                   /* ADDED BY DPD */
  sma->tlci_flows[6] = bfcc;                                  /* ADDED BY DPD */

  /*    COMPUTE TOTAL CHANNEL INFLOW FOR THE TIME INTERVAL.  */
  sma->tlci = roimp + sdro + ssur + sif + bfcc;

  /*    COMPUTE E4-ET FROM RIPARIAN VEGETATION. */
  e4 = (edmnd - eused) * sma->riva;

  /*    SUBTRACT E4 FROM CHANNEL INFLOW */
  sma->tlci -= e4;

  if (sma->tlci < 0.) {
    e4 += sma->tlci;
    sma->tlci = 0.;
  }
  fsum1->srot += sma->tlci;

  /*    COMPUTE TOTAL EVAPOTRANSPIRATION-TET */
  eused *= parea;
  tet = eused + e5 + e4;
  fsum1->sett += tet;
  fsum1->se1 += (e1 * parea);
  fsum1->se3 += (e3 * parea);
  fsum1->se4 += e4;
  fsum1->se5 += e5;

  /*   CHECK THAT ADIMC >= UZTWC */
  if (sma->adimc < sma->uztwc)
    sma->adimc = sma->uztwc;

  /*   ADD TO SUMS OF RUNOFF COMPONENTS. */
  rsum[0] += sma->tlci;
  rsum[1] += roimp;
  rsum[2] += sdro;
  rsum[3] += ssur;
  rsum[4] += sif;
  rsum[5] += bfs;
  rsum[6] += bfp;

  /*      if (iter1>=5050) printf("end of fland1\n");*/
  /*      printf("end of fland\n");*/
  /*   END OF FUNCTION FLAND1 */
}

#undef EPS
