## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

#' Typical initial model used in Data-Based Mechanistic modelling.
#'
#' Typical initial model used in Data-Based Mechanistic modelling.  Rainfall is
#' scaled by corresponding streamflow values raised to a power.  This SMA uses
#' streamflow data, so can not be used for prediction.
#'
#' @name dbm
#' @aliases dbm dbm.sim absorbScale.hydromad.dbm
#' @param DATA time-series-like object with columns \code{P} (precipitation)
#' and \code{Q} (streamflow).
#' @param power power to apply to streamflow values.
#' @param qlag number of time steps to lag the streamflow (relative to
#' rainfall) before multiplication.
#' @param scale constant multiplier of the result, for mass balance.  If this
#' parameter is set to \code{NA} (as it is by default) in
#' \code{\link{hydromad}} it will be set by mass balance calculation.
#' @param return_state ignored.
#' @return the simulated effective rainfall, a time series of the same length
#' as the input series.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "dbm")} to work with models as objects
#' (recommended).
#' @keywords models
#' @examples
#'
#' ## view default parameter ranges:
#' str(hydromad.options("dbm"))
#'
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "dbm", routing = "expuh")
#' mod0
#'
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, power = 0.5, qlag = 0, tau_s = 10)
#'
#' xyplot(cbind(HydroTestData, dbm.Q = predict(mod1)))
#'
#' ## show effect of increase/decrease in each parameter
#' parRanges <- list(power = c(0.01, 0.9), qlag = c(-1, 2))
#' parsims <- mapply(
#'   val = parRanges, nm = names(parRanges),
#'   FUN = function(val, nm) {
#'     lopar <- min(val)
#'     hipar <- max(val)
#'     names(lopar) <- names(hipar) <- nm
#'     fitted(runlist(
#'       decrease = update(mod1, newpars = lopar),
#'       increase = update(mod1, newpars = hipar)
#'     ))
#'   }, SIMPLIFY = FALSE
#' )
#'
#' xyplot.list(parsims,
#'   superpose = TRUE, layout = c(1, NA),
#'   main = "Simple parameter perturbation example"
#' ) +
#'   latticeExtra::layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
#' @export
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


#' @export
absorbScale.hydromad.dbm <- function(object, gain, ...) {
  absorbScale.hydromad.scalar(object, gain, parname = "scale")
}
