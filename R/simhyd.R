## hydromad: Hydrological Modelling and Analysis of Data
## Rewrite
## willem 20151014
## Coded based on diagram and description in:
# Chiew et al 2009 WATER RESOURCES RESEARCH, VOL. 45, W10414, doi:10.1029/2008WR007338, 2009

## SimHyd model
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

# Routing based on Muskinghum
simhydrouting.sim <- function(U, DELAY = 1, X_m = 0.2,
                              epsilon = hydromad.getOption("sim.epsilon"),
                              return_components = FALSE) {
  X <- rep(0, length(U))
  inAttr <- attributes(U)
  U <- as.ts(U)
  bad <- is.na(U)
  U[bad] <- 0
  if (2 * DELAY * X_m < 1 & 2 * DELAY * (1 - X_m) > 1) {
    # Muskingum components
    C0 <- (-DELAY * X_m + 0.5) / (DELAY * (1 - X_m) + 0.5)
    # print(C0)
    C1 <- (DELAY * X_m + 0.5) / (DELAY * (1 - X_m) + 0.5)
    # print(C1)
    C2 <- (DELAY * (1 - X_m) - 0.5) / (DELAY * (1 - X_m) + 0.5)
    # print(C2)
  } else {
    C0 <- 0
    C1 <- 1
    C2 <- 0
    # print("model parameters adjusted")
  }

  if (C0 + C1 + C2 != 1) {
    C0 <- 0
    C1 <- 1
    C2 <- 0
    # print("model parameters adjusted again")
  }
  # print(C0+C1+C2)
  # if (round(C0+C1+C2)!=1)  C0 <- 0; C1 <- 1; C2 <- 0

  X[1] <- U[1]
  for (t in 1:(length(U) - 1)) {
    X[t + 1] <- C0 * U[t + 1] + C1 * U[t] + C2 * X[t]
    # print(X[t+1])
  }
  X[abs(X) < epsilon] <- 0
  X[bad] <- NA
  attributes(X) <- inAttr
  X
}

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

simhydrouting.ranges <- function() {
  list(
    DELAY = c(0.1, 5),
    X_m = c(0.01, 0.5)
  )
}
