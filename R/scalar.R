## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


#' Simple constant runoff proportion
#'
#' Simple constant runoff proportion: a constant fraction of rainfall reaches
#' the stream.
#'
#' @name scalar
#' @aliases scalar.sim absorbScale.hydromad.scalar
#' @param DATA time-series-like object with a column P (precipitation).
#' @param scale fraction of rainfall that becomes effective.  If this parameter
#' is set to \code{NA} (as it is by default) in \code{\link{hydromad}} it will
#' be set by mass balance calculation.
#' @param return_state ignored.
#' @return the simulated effective rainfall, a time series of the same length
#' as the input series.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "scalar")} to work with models as
#' objects (recommended).
#' @keywords models
#' @examples
#'
#' ## view default parameter ranges:
#' str(hydromad.options("scalar"))
#'
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "scalar", routing = "expuh")
#' mod0
#'
#' ## simulate with some arbitrary parameter values
#' testQ <- predict(update(mod0, scale = 0.5, tau_s = 10))
#' xyplot(cbind(HydroTestData[, 1:2], scalar.Q = testQ))
#' @export
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

#' @importFrom stats coef
#' @export
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
