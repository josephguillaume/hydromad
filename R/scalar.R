## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

scalar.sim <-
  function(DATA, scale, return_state = FALSE) {
    if (NCOL(DATA) > 1) stopifnot("P" %in% colnames(DATA))
    P <- if (NCOL(DATA) > 1) DATA[, "P"] else DATA
    ## special value scale = NA used for initial run
    if (is.na(scale)) {
      scale <- 1
    }
    ## check values
    stopifnot(scale >= 0)
    ## compute effective rainfall U
    scale * P
  }

scalar.ranges <- function() {
  list(scale = NA_real_)
}

absorbScale.hydromad.scalar <- function(object, gain, parname = "scale", ...) {
  if (gain <= 0) {
    return(NULL)
  }
  coeff <- coef(object, which = "sma")
  if (parname %in% names(coeff) == FALSE) {
    return(NULL)
  }
  scale <- coeff[[parname]]
  ## we only want to do this when scale is NA (special value)
  if (is.null(scale) || !is.na(scale)) {
    return(NULL)
  }
  scale <- 1
  scale <- scale * gain
  object$parlist[[parname]] <- scale
  object$call[[parname]] <- signif(scale, 6)
  object$U <- scale * object$U
  object
}
