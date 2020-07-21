## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


#' Simple time-varying runoff proportion
#'
#' Simple time-varying runoff proportion. Rainfall is scaled by the runoff
#' coefficient estimated in a moving window.  This SMA uses streamflow data, so
#' can not be used for prediction.
#'
#'
#' @name runoffratio
#' @aliases runoffratio.sim absorbScale.hydromad.runoffratio
#' @param DATA time-series-like object with columns \code{P} (precipitation)
#' and \code{Q} (streamflow).
#' @param width width of the time window (in time steps) in which to estimate
#' the runoff coefficient.
#' @param kernel type of window used to estimate the runoff coefficient: 1 is
#' rectangular, 2 is triangular-weighted, 3 is Gaussian-like.
#' @param sides 2 for time-centered estimates, 1 for estimates using data
#' backward in time only.
#' @param rrthresh a theshold value of the runoff ratio, below which there is
#' no effective rainfall.
#' @param qlag number of time steps to lag the streamflow (relative to
#' rainfall) before estimating the runoff coefficient.
#' @param scale constant multiplier of the result, for mass balance.  If this
#' parameter is set to \code{NA} (as it is by default) in
#' \code{\link{hydromad}} it will be set by mass balance calculation.
#' @param return_state ignored.
#' @return the simulated effective rainfall, a time series of the same length
#' as the input series.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "runoffratio")} to work with models as
#' objects (recommended).
#' @keywords models
#' @examples
#'
#' ## view default parameter ranges:
#' str(hydromad.options("runoffratio"))
#'
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "runoffratio", routing = "expuh")
#' mod0
#'
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, width = 30, rrthresh = 0.2, tau_s = 10)
#'
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[, 1:2], runoffratio = testQ))
#'
#' ## show effect of increase/decrease in each parameter
#' parRanges <- list(
#'   width = c(10, 180), qlag = c(-30, 30),
#'   rrthresh = c(0, 0.5)
#' )
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
runoffratio.sim <-
  function(DATA, width = 30, kernel = 2, sides = 2, rrthresh = 0,
           qlag = 0, scale = 1, return_state = FALSE) {
    stopifnot(c("P", "Q") %in% colnames(DATA))
    P <- DATA[, "P"]
    Q <- DATA[, "Q"]
    ## special value scale = NA used for initial run for scaling
    if (is.na(scale)) {
      scale <- 1
    }
    ## check values
    stopifnot(kernel >= 1)
    ## synchronise Q lagged by qlag
    Q <- shiftWindow(Q, round(qlag), and.lag = TRUE)
    ## compute effective rainfall U
    ## estimate U as scaled P, scaled in a moving window
    ## (runoff coefficient)
    width <- round(width)
    sm <- simpleSmoothTs(cbind(Q, P), width = width, c = kernel, sides = sides)
    rr <- sm[, "Q"] / sm[, "P"]
    rr <- na.locf(na.locf(rr, na.rm = FALSE), fromLast = TRUE, na.rm = FALSE)
    rr[!is.finite(rr)] <- 0
    rr[rr < rrthresh] <- 0
    U <- scale * P * rr

    ans <- U
    if (return_state) {
      ans <- cbind(U = U, rr = rr)
    }
    return(ans)
  }



runoffratio.ranges <- function() {
  list(
    rrthresh = c(0, 0.2),
    scale = NA
  )
}


#' @export
absorbScale.hydromad.runoffratio <-
  function(object, gain, ...) {
    absorbScale.hydromad.scalar(object, gain, parname = "scale")
  }
