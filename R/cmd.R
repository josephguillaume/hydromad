## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

#' IHACRES Catchment Moisture Deficit (CMD) model
#'
#' The Catchment Moisture Deficit (CMD) effective rainfall model for IHACRES.
#' It is a conceptual-type model, where input rainfall is partitioned
#' explicitly into drainage, evapo-transpiration, and changes in catchment
#' moisture.
#'
#' The mass balance step is: \deqn{M[t] = M[t-1] - P[t] + E_T[t] + U[t]}
#'
#' where \eqn{M} represents catchment moisture deficit (CMD), constrained below
#' by 0 (the nominal fully saturated level).  P is catchment areal rainfall,
#' \eqn{E_T} is evapo-transpiration, and U is drainage (effective rainfall).
#' All are, typically, in units of mm per time step.
#'
#' Rainfall effectiveness (i.e. drainage proportion) is a simple
#' \emph{instantaneous} function of the CMD, with a threshold at \eqn{M = d}.
#' In the default linear form this is:
#'
#' \deqn{\frac{\mathrm{d}U}{\mathrm{d}P} = 1 - \min(1, M/d)}{ dU/dP = 1 -
#' min(1, M/d)}
#'
#' The trigonometric form is
#'
#' \deqn{\frac{\mathrm{d}U}{\mathrm{d}P} = 1 - \min(1, \sin^2(\pi M / 2d))}{
#' dU/dP = 1 - min(1, sin^2(pi M / 2d))}
#'
#' The power form is
#'
#' \deqn{\frac{\mathrm{d}U}{\mathrm{d}P} = 1 - \min(1, (M/d)^a)}{ dU/dP = 1 -
#' min(1, (M/d)^a)} where a = 10 ^ (shape / 50)
#'
#' The actual drainage each time step involves the integral of these relations.
#'
#' Evapo-transpiration is also a simple function of the CMD, with a threshold
#' at \eqn{M = f d}{M = f * d}: \deqn{E_T[t] = e E[t] \min(1,
#' \exp\left(2\left(1 - \frac{M_f}{fd}\right)\right))}{ E_T[t] = e E[t] \min(1,
#' \exp(2(1 - M_f / (fd))))}
#'
#' Note that the evapo-transpiration calculation is based on \eqn{M_f}, which
#' is the CMD after precipitation and drainage have been accounted for.
#'
#' @name IHACRES.CMD.model
#' @aliases cmd cmd.sim
#' @param DATA a \code{\link{ts}}-like object with named columns: \describe{
#' \item{list("P")}{ time series of areal rainfall depths, usually in mm. }
#' \item{list("E")}{ time series of potential evapo-transpiration, or more
#' typically, temperature as an indicator of this. } }
#' @param f CMD stress threshold as a proportion of \code{d}.
#' @param e temperature to PET conversion factor.
#' @param d CMD threshold for producing flow.
#' @param shape defines form of the \eqn{dU/dP} relationship: \code{shape = 0}
#' is the linear form, \code{shape = 1} is the trigonometric form, and
#' \code{shape > 1} is the power form.
#' @param M_0 starting CMD value.
#' @param return_state to return state variables as well as the effective
#' rainfall.
#' @return \code{cmd.sim} returns the modelled time series of effective
#' rainfall, or if \code{return_state = TRUE}, a multi-variate time series with
#' named columns \code{U} (effective rainfall), \code{CMD} and \code{ET}
#' (evapo-transpiration \eqn{E_T}).
#' @note Normally compiled C code is used for simulation, but if
#' \code{return_state = TRUE} a slower implementation in R is used.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "cmd")} to work with models as objects
#' (recommended).
#' @references Croke, B.F.W. and A.J. Jakeman (2004), A Catchment Moisture
#' Deficit module for the IHACRES rainfall-runoff model, \emph{Environmental
#' Modelling and Software}, 19(1): 1-5.
#'
#' Croke, B.F.W. and A.J. Jakeman (2005), Corrigendum to ``A Catchment Moisture
#' Deficit module for the IHACRES rainfall-runoff model'' [Environ. Model.
#' Softw. 19 (1) (2004) 1-5], \emph{Environmental Modelling and Software},
#' 20(7): 977.
#' @keywords models
#' @examples
#'
#' ## view default parameter ranges:
#' str(hydromad.options("cmd"))
#'
#' data(Canning)
#' x <- cmd.sim(Canning[1:1000, ],
#'   d = 200, f = 0.7, e = 0.166,
#'   return_state = TRUE
#' )
#' xyplot(x)
#'
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "cmd", routing = "expuh")
#' mod0
#'
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, d = 200, f = 0.5, e = 0.1, tau_s = 10)
#'
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[, 1:2], cmd = testQ))
#'
#' ## show effect of increase/decrease in each parameter
#' parlist <- list(
#'   d = c(50, 550), f = c(0.01, 3),
#'   e = c(0.01, 1.5)
#' )
#' parsims <- mapply(
#'   val = parlist, nm = names(parlist),
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
#' @rdname IHACRES.CMD.model
#' @useDynLib hydromad sma_cmd
#' @export
cmd.sim <-
  function(DATA,
           f, e, d, shape = 0,
           M_0 = d / 2,
           return_state = FALSE) {
    stopifnot(c("P", "E") %in% colnames(DATA))
    ## check values
    stopifnot(f >= 0)
    stopifnot(e >= 0)
    stopifnot(d >= 0)
    stopifnot(shape >= 0)
    ## default initial state
    if (is.na(M_0)) M_0 <- d / 2

    ## f is expressed as a proportion of d
    g <- f * d

    inAttr <- attributes(DATA[, 1])
    DATA <- as.ts(DATA)
    P <- DATA[, "P"]
    E <- DATA[, "E"]
    ## skip over missing values (maintaining the state M)
    bad <- is.na(P) | is.na(E)
    P[bad] <- 0
    E[bad] <- 0
    ## TODO: return state from C code
    COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
    if (COMPILED) {
      ans <- .C(sma_cmd,
        as.double(P),
        as.double(E),
        as.integer(NROW(DATA)),
        as.double(g),
        as.double(e),
        as.double(d),
        as.double(shape),
        as.double(M_0),
        U = double(NROW(DATA)),
        M = double(NROW(DATA)),
        ET = double(NROW(DATA)),
        NAOK = FALSE, PACKAGE = "hydromad"
      )
      U <- ans$U
      M <- ans$M
      ET <- ans$ET
    } else {
      ## implementation in R for cross-checking (slow)
      U <- M <- ET <- P
      M_prev <- M_0
      for (t in seq(1, length(P))) {
        ## default, for when P[t] == 0:
        Mf <- M_prev
        ## select form of dU/dP relationship
        if (P[t] > 0) {
          ## rainfall reduces CMD (Mf)
          if (shape < 1) {
            ## linear form: dU/dP = 1 - (M/d) for M < d
            if (M_prev < d) {
              Mf <- M_prev * exp(-P[t] / d)
            } else if (M_prev < d + P[t]) {
              Mf <- d * exp((-P[t] + M_prev - d) / d)
            } else {
              Mf <- M_prev - P[t]
            }
          }
          else if (shape == 1) {
            ## trigonometric form: dU/dP = 1 - sin^2(pi M / 2d) for M < d
            if (M_prev < d) {
              Mf <- 1 / tan((M_prev / d) * (pi / 2))
              Mf <- (2 * d / pi) * atan(1 / (pi * P[t] / (2 * d) + Mf))
            } else if (M_prev < d + P[t]) {
              Mf <- (2 * d / pi) * atan(2 * d / (pi * (d - M_prev + P[t])))
            } else {
              Mf <- M_prev - P[t]
            }
          }
          else { ## shape > 1
            ## power form: dU/dP = 1 - (M/d)^a for M < d
            a <- 10^(shape / 50)
            if (M_prev < d) {
              Mf <- M_prev * (1 - ((1 - a) * P[t] / (d^a)) /
                (M_prev^(1 - a)))^(1 / (1 - a))
            } else if (M_prev < d + P[t]) {
              Mf <- d * (1 - (1 - a) * (P[t] - M_prev + d) / d)^(1 / (1 - a))
            } else {
              Mf <- M_prev - P[t]
            }
          }
        }
        ## drainage (rainfall not accounted for in -dM)
        U[t] <- max(0, P[t] - M_prev + Mf)
        ## evapo-transpiration
        ET[t] <- e * E[t] * min(1, exp(2 * (1 - Mf / g)))
        ET[t] <- max(0, ET[t])
        ## mass balance
        M[t] <- M_prev - P[t] + U[t] + ET[t]
        M_prev <- M[t] <- max(0, M[t])
      }
    }
    attributes(U) <- inAttr
    ## re-insert missing values
    U[bad] <- NA
    ans <- U
    if (return_state) {
      attributes(M) <- attributes(ET) <- attributes(U)
      ans <- cbind(U = U, CMD = M, ET = ET)
    }
    return(ans)
  }

cmd.ranges <- function() {
  list(
    f = c(0.01, 3),
    e = c(0.01, 1.5),
    d = c(50, 550),
    shape = 0
  )
}
