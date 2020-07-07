## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

dbm.sim <-
  function(DATA, power, qlag = 0, scale = 1, return_state = FALSE) {
    stopifnot(c("P", "Q") %in% colnames(DATA))
    P <- DATA[, "P"]
    Q <- DATA[, "Q"]
    ## special value scale = NA used for initial run for scaling
    if (is.na(scale)) {
      scale <- 1
    }
    ## check values
    stopifnot(scale >= 0)
    ## synchronise Q lagged by qlag
    Qd <- shiftWindow(Q, round(qlag), and.lag = TRUE)
    ## compute effective rainfall U
    scale * P * Qd^power
  }

dbm.ranges <- function() {
  list(
    power = c(0, 0.9),
    qlag = c(-1, 2),
    scale = NA_real_
  )
}

absorbScale.hydromad.dbm <- function(object, gain, ...) {
  absorbScale.hydromad.scalar(object, gain, parname = "scale")
}
