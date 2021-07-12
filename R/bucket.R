## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## Bucket-type Soil Moisture Accounting models.
## Bai et. al. (2009), Environmental Modelling and Software.
## Model S2.

#' Single-bucket Soil Moisture Accounting models
#'
#' Single-bucket Soil Moisture Accounting models with saturated/unsaturated
#' zones and interception.
#'
#' From formulations given in Bai et. al. (2009), which were based on Farmer
#' et. al. (2003).
#'
#' The general mass balance structure is: \deqn{dS/dt = p - q(S) - e(S, Ep)}
#'
#' The default parameter ranges were also taken from Bai et. al. (2009).
#'
#' @name bucket
#' @aliases bucket bucket.sim
#' @param DATA time-series-like object with columns P (precipitation, mm) and E
#' (potential evapo-transpiration, mm).
#' @param Sb Maximum soil water storage (mm).
#' @param fc Field capacity (0 - 1).
#' @param a.ei Interception coefficient (\eqn{\alpha_{ei}}).
#' @param M Fraction of catchment area covered by deep rooted vegetation.
#' @param a.ss Recession coefficients for subsurface flow from saturated zone
#' (\eqn{\alpha_{ss}}).
#' @param etmult Multiplier for the \code{E} input data.
#' @param S_0 Initial soil moisture level as a fraction of \code{Sb}.
#' @param return_state to return the series U, S and ET (evapotranspiration).
#' @return the simulated effective rainfall, a time series of the same length
#' as the input series.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "bucket")} to work with models as
#' objects (recommended).
#' @references Farmer, D., M. Sivapalan, Farmer, D. (2003). Climate, soil and
#' vegetation controls upon the variability of water balance in temperate and
#' semiarid landscapes: downward approach to water balance analysis.
#' \emph{Water Resources Research} 39(2), p 1035.
#'
#' Bai, Y., T. Wagener, P. Reed (2009). A top-down framework for watershed
#' model evaluation and selection under uncertainty. \emph{Environmental
#' Modelling and Software} 24(8), pp. 901-916.
#' @keywords models
#' @examples
#'
#' ## view default parameter ranges:
#' str(hydromad.options("bucket"))
#'
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "bucket", routing = "expuh")
#' mod0
#'
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0,
#'   Sb = 10, fc = 0.5, M = 0.5, etmult = 0.05,
#'   a.ei = 0.05, a.ss = 0.01, tau_s = 10
#' )
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[, 1:2], bucket = testQ))
#'
#' ## show effect of increase/decrease in each parameter
#' parRanges <- hydromad.getOption("bucket")
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
#'   strip = FALSE, strip.left = TRUE,
#'   main = "Simple parameter perturbation example"
#' ) +
#'   latticeExtra::layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
#' @useDynLib hydromad sma_bucket
#' @export
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
        NAOK = FALSE, PACKAGE = "hydromad"
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
        ET[t] <- Eintc + min(S[t], Etrans + Ebare)
        ## mass balance
        S[t] <- S_prev + P[t] - ET[t]
        ## drainage (saturation excess)
        Use <- max(0, S[t] - Sb)
        S[t] <- S[t] - Use
        ## drainage (sub-surface)
        Uss <- max(0, a.ss * (S[t] - Sfc))
        S[t] <- S[t] - Uss
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
