## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


## Australian Water Balance Model.
## Boughton (2004), Environmental Modelling and Software.
## This is the bucket component, to be used with
## a two-store routing function:
## routing = "expuh", rfit = list("sriv", order = c(2,1))
awbm.sim <-
  function(DATA, cap.ave,
           cap1 = 0.01 * cap.ave / area1,
           cap2 = 0.33 * cap.ave / area2,
           cap3 = 0.66 * cap.ave / area3,
           area2 = 0.433, area3 = 0.433,
           area1 = 1 - area2 - area3,
           etmult = 1,
           S1_0 = 0, S2_0 = 0, S3_0 = 0,
           return_state = FALSE) {
    stopifnot(c("P", "E") %in% colnames(DATA))
    ## check values
    stopifnot(cap1 >= 0)
    stopifnot(cap2 >= 0)
    stopifnot(cap3 >= 0)
    stopifnot(S1_0 >= 0)
    stopifnot(S2_0 >= 0)
    stopifnot(S3_0 >= 0)
    stopifnot(0 <= area1)
    stopifnot(0 <= area2)
    stopifnot(0 <= area3)
    stopifnot(area1 == 1 - area2 - area3)

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
      ans <- .C(sma_awbm,
        as.double(P),
        as.double(E),
        as.integer(NROW(DATA)),
        as.double(cap1),
        as.double(cap2),
        as.double(cap3),
        as.double(area1),
        as.double(area2),
        as.double(area3),
        as.double(S1_0),
        as.double(S2_0),
        as.double(S3_0),
        U = double(NROW(DATA)),
        S1 = double(NROW(DATA)),
        S2 = double(NROW(DATA)),
        S3 = double(NROW(DATA)),
        NAOK = FALSE, DUP = FALSE, PACKAGE = "hydromad"
      )
      U <- ans$U
      S1 <- ans$S1
      S2 <- ans$S2
      S3 <- ans$S3
    } else {
      ## implementation in R for cross-checking (slow)
      U <- S1 <- S2 <- S3 <- P
      S1_prev <- S1_0
      S2_prev <- S2_0
      S3_prev <- S3_0
      for (t in seq(1, length(P))) {
        ## rainfall
        S1[t] <- S1_prev + P[t]
        S2[t] <- S2_prev + P[t]
        S3[t] <- S3_prev + P[t]
        ## saturation excess
        U1 <- max(0, S1[t] - cap1)
        U2 <- max(0, S2[t] - cap2)
        U3 <- max(0, S3[t] - cap3)
        U[t] <- area1 * U1 + area2 * U2 + area3 * U3
        ## evapotranspiration
        S1_prev <- S1[t] <- max(0, S1[t] - U1 - E[t])
        S2_prev <- S2[t] <- max(0, S2[t] - U2 - E[t])
        S3_prev <- S3[t] <- max(0, S3[t] - U3 - E[t])
      }
    }
    attributes(U) <- inAttr
    ## re-insert missing values
    U[bad] <- NA
    ans <- U
    if (return_state) {
      attributes(S1) <- attributes(S2) <- attributes(S3) <- attributes(U)
      ans <- cbind(U = U, S1 = S1, S2 = S2, S3 = S3)
    }
    return(ans)
  }

awbm.ranges <- function() {
  list(
    cap.ave = c(1, 1000),
    etmult = c(0.01, 1)
  )
}
