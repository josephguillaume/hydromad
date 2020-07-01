## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## Bucket-type Soil Moisture Accounting models.
## Bai et. al. (2009), Environmental Modelling and Software.
## Model S2.
bucket.sim <-
  function(DATA,
           Sb, fc = 1, a.ei = 0, M = 0, a.ss = 0,
           etmult = 1, S_0 = 0.5,
           return_state = FALSE) {
    stopifnot(c("P", "E") %in% colnames(DATA))
    ## check values
    stopifnot(Sb > 0)
    stopifnot(S_0 >= 0)
    stopifnot(0 < fc && fc <= 1)
    stopifnot(0 <= a.ei && a.ei <= 1)
    stopifnot(0 <= M && M <= 1)
    stopifnot(0 <= a.ss && a.ss <= 1)

    ## fc is expressed as a proportion of Sb
    Sfc <- fc * Sb
    ## as is S_0
    S_0 <- S_0 * Sb

    inAttr <- attributes(DATA[, 1])
    DATA <- as.ts(DATA)
    P <- DATA[, "P"]
    E <- DATA[, "E"] * etmult

    ## skip over missing values (maintaining the state S)
    bad <- is.na(P) | is.na(E)
    P[bad] <- 0
    E[bad] <- 0
    COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
    if (COMPILED) {
      ans <- .C(sma_bucket,
        as.double(P),
        as.double(E),
        as.integer(NROW(DATA)),
        as.double(Sb),
        as.double(fc),
        as.double(Sfc),
        as.double(a.ei),
        as.double(M),
        as.double(a.ss),
        as.double(S_0),
        U = double(NROW(DATA)),
        S = double(NROW(DATA)),
        ET = double(NROW(DATA)),
        NAOK = FALSE, DUP = FALSE, PACKAGE = "hydromad"
      )
      U <- ans$U
      S <- ans$S
      ET <- ans$ET
    } else {
      ## implementation in R for cross-checking (slow)
      U <- S <- ET <- P
      S_prev <- S_0
      for (t in seq(1, length(P))) {
        ## evapo-transpiration
        Eintc <- a.ei * P[t]
        S[t] <- min(Sb, S_prev + P[t] - Eintc)
        Etrans <- M * min(1, S[t] / Sfc) * E[t]
        Ebare <- (1 - M) * (S[t] / Sb) * E[t]
        ET[t] <- Eintc + Etrans + Ebare
        ## mass balance
        S[t] <- max(0, S_prev + P[t] - ET[t])
        ## drainage (saturation excess)
        Use <- max(0, S[t] - Sb)
        S[t] <- max(0, S[t] - Use)
        ## drainage (sub-surface)
        Uss <- max(0, a.ss * (S[t] - Sfc))
        S[t] <- max(0, S[t] - Uss)
        U[t] <- Use + Uss
        S_prev <- S[t]
      }
    }
    attributes(U) <- inAttr
    ## re-insert missing values
    U[bad] <- NA
    ans <- U
    if (return_state) {
      attributes(S) <- attributes(ET) <- attributes(U)
      ans <- cbind(U = U, S = S, ET = ET)
    }
    return(ans)
  }

bucket.ranges <- function() {
  list(
    Sb = c(0.1, 1200),
    fc = c(0.01, 1),
    a.ei = c(0, 0.49),
    M = c(0, 1),
    a.ss = c(0, 0.5)
  )
}
