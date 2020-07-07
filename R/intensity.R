## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

intensity.sim <-
  function(DATA, power, maxP = 500, scale = 1, return_state = FALSE) {
    if (NCOL(DATA) > 1) stopifnot("P" %in% colnames(DATA))
    P <- if (NCOL(DATA) > 1) DATA[, "P"] else DATA
    ## special value scale = NA used for initial run for scaling
    if (is.na(scale)) {
      scale <- 1
    }
    ## check values
    stopifnot(power >= 0)
    stopifnot(maxP > 0)
    stopifnot(scale >= 0)
    ## compute effective rainfall U
    scale * P * pmin((P^power) / (maxP^power), 1)
  }

intensity.ranges <- function() {
  list(
    power = c(0, 2),
    maxP = c(100, 1000),
    scale = NA_real_
  )
}

absorbScale.hydromad.intensity <- function(object, gain, ...) {
  absorbScale.hydromad.scalar(object, gain, parname = "scale")
}
