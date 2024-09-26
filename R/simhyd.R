## hydromad: Hydrological Modelling and Analysis of Data
## SimHyd
## Coded based on diagram and description in:
# Chiew et al 2009 WATER RESOURCES RESEARCH, VOL. 45, W10414, doi:10.1029/2008WR007338, 2009

#' SimHyd model
#'
#' Coded based on diagram and description in:
#' Chiew et al 2009 WATER RESOURCES RESEARCH, 
#' VOL. 45, W10414, doi:10.1029/2008WR007338, 2009
#'
#' @name simhyd
#' @aliases simhyd.sim simhyd.ranges
#'
#' @param DATA time-series-like object with columns P (precipitation, mm) and E
#' (potential evapo-transpiration, mm).
#' @param U effective rainfall series.
#' @param INSC Interception storage capacity (mm)
#' @param COEFF maximum infiltration loss (mm)
#' @param SQ Infiltration loss exponent
#' @param SMSC Soil moisture store capacity (mm)
#' @param SUB Constant of proportionality in interflow equation
#' @param CRAK Constant of proportionality in groundwater recharge equation
#' @param K Baseflow linear recession parameter
#' @param etmult mutliplier to convert daily maximum temperature to estimates
#'  of potential evaporation
#' @param GWt0 Initial state of the groundwater store (mm)
#' @param SMSt0 Initial state of soil moisture store (mm)
#' @param return_state to return the series U, S (storage) and ET
#' (evapotranspiration).
#' @author Willem Vervoort \email{willemvervoort@@gmail.com} and Joseph Guillaume
#' \email{josephguillaume@@gmail.com}
#' @seealso \code{\link{hydromad}(sma = "simhyd", routing = "simhydrouting")} to
#' work with models as objects (recommended).
#' @references Chiew, F. H. S., Teng, J., Vaze, J., Post, D. A., Perraud, J. M.,
#' Kirono, D. G. C., and Viney, N. R. (2009), 
#' Estimating climate change impact on runoff across southeast Australia: 
#' Method, results, and implications of the modeling method, 
#' \emph{Water Resour. Res.}, 45, W10414, \url{doi:10.1029/2008WR007338}.
#'
#' @keywords models
#' @examples
#'
#' ## view default parameter ranges:
#' str(hydromad.getOption("simhyd"))
#'
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "simhyd", routing = "simhydrouting")
#' mod0
#'
#' data(Cotter)
#' x <- Cotter[1:1000]
#'
#' # Specify simhyd model
#' mod0 <- hydromad(x, sma = "simhyd", routing = "simhydrouting")
#' # specify parameter ranges for a few parameters
#' mod0 <- update(mod0,  COEFF=c(0,400), SQ=c(0.1,5), 
#' K=c(0,1))
#' # Allow etmult to vary, because we're using temperature data instead of PET.
#' mod0 <- update(mod0, etmult = c(0.05, 1.5))
#' 
#' mod0
#'
#' ## now try to fit the free parameters
#' set.seed(10)
#' fit1 <- fitByOptim(mod0)
#'
#' fit1
#' summary(fit1)
#' xyplot(fit1)


#' @useDynLib hydromad _hydromad_simhyd_sim
#' @useDynLib hydromad
#' @importFrom Rcpp sourceCpp
#' @export
simhyd.sim <-
  function(DATA,
           INSC, COEFF,
           SQ,
           SMSC, SUB, CRAK, K,
           etmult = 0.15,
           GWt0 = 0, SMSt0 = 0.5,
           return_state = FALSE)

           # See Figure 2 in Chiew et al. 2009
           # INSC interception store capacity (mm)
           # COEFF maximum infiltration loss
           # SQ Infiltration loss exponent
           #  SMSC = Soil Moisture Storage Capacity
           # SUB constant of proportionality in interflow equation
           # CRAK constant of proportionality in groundwater rechareg equation
           # K baseflow linear recession parameter
  # etmult = added parameter to convert maxT to PET

  {
    stopifnot(c("P", "E") %in% colnames(DATA))
    ## check values
    stopifnot(INSC >= 0)
    stopifnot(COEFF >= 0)
    stopifnot(SQ >= 0)
    stopifnot(SMSC >= 0)
    stopifnot(SUB >= 0)
    stopifnot(CRAK >= 0)
    stopifnot(K >= 0)
    xpar <-
      c(INSC, COEFF, SQ, SMSC, SUB, CRAK, K)

    inAttr <- attributes(DATA[, 1])
    DATA <- as.ts(DATA)

    P <- DATA[, "P"]
    E <- etmult * DATA[, "E"]
    ## skip over missing values
    bad <- is.na(P) | is.na(E)
    P[bad] <- 0
    E[bad] <- 0

    COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
    if (COMPILED) {
      # run the cpp version
      ans <- simhyd_sim(
        P, E, INSC, COEFF, SQ, SMSC,
        SUB, CRAK, K, GWt0, SMSt0
      )
      U <- ans$U
      ET <- ans$ET
    } else { ## very slow, even on my x64
      U <- IMAX <- INT <- INR <- RMO <- IRUN <-
        ET <- SRUN <- REC <- SMF <- POT <- BAS <-
        SMS <- GW <- rep(NA_real_, length(P))
      GWt1 <- GWt0
      SMSt1 <- SMSt0 * SMSC
      # run through a loop
      for (t in seq(1, length(P))) {
        # interception store
        IMAX[t] <- min(INSC, E[t])
        # print(IMAX[t])
        # calculate interception
        INT[t] <- min(IMAX[t], P[t])
        # print(INT[t])
        # calculate interception runoff (INR)
        INR[t] <- P[t] - INT[t]
        # print(INR[t])
        # Calculate infiltration capacity
        RMO[t] <- min(COEFF * exp(-SQ * SMSt1 / SMSC), INR[t])
        # print(RMO[t])
        # calculate direct runoff
        IRUN[t] <- INR[t] - RMO[t]
        # print(IRUN[t])
        # SRUN (Saturation excess runoff and interflow)
        SRUN[t] <- SUB * SMSt1 / SMSC * RMO[t]
        # print(SRUN[t])
        # calculate Recharge
        REC[t] <- CRAK * SMSt1 / SMSC * (RMO[t] - SRUN[t])
        # print(REC[t])
        # INfiltration into soil store (SMF)
        SMF[t] <- RMO[t] - SRUN[t] - REC[t]
        # print(SMF[t])
        # calculate potential ET
        POT[t] <- E[t] - INT[t]
        # Calculate Soil ET
        ET[t] <- min(10 * SMSt1 / SMSC, POT[t])
        # print(ET[t])
        # calculate SMS overflow (see Figure 2 in Chiew et al 2009)
        # calculate soil moisture storage
        SMS[t] <- SMSt1 + SMF[t] - ET[t]
        if (SMS[t] > SMSC) {
          SMS[t] <- SMSC
          REC[t] <- REC[t] + SMS[t] - SMSC
        }
        SMSt1 <- SMS[t]
        # calculate baseflow
        BAS[t] <- K * GWt1
        # Calculate GW storage
        GW[t] <- GWt1 + REC[t] - BAS[t]
        GWt1 <- GW[t]
        # Calculate runoff
        U[t] <- IRUN[t] + SRUN[t] + BAS[t]
      }
    }
    ## make it a time series object again
    attributes(U) <- inAttr
    attributes(ET) <- inAttr
    ## re-insert missing values
    U[bad] <- NA
    ET[bad] <- NA
    if (return_state == TRUE) {
      return(merge(U = U, ET = ET))
    } else {
      return(U)
    }
  }


#' @export
simhyd.ranges <- function() {
  list(
    INSC = c(0, 50),
    COEFF = c(0.0, 400),
    SQ = c(0, 10),
    SMSC = c(1, 1000),
    SUB = c(0.0, 1),
    CRAK = c(0.0, 1),
    K = c(0.0, 1),
    etmult = c(0.01, 1)
  )
}

