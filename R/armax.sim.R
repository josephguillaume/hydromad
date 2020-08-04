## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

#' ARMAX Transfer Function models
#'
#' ARMAX linear transfer functions with a single input and single output
#' series. Can be used as a general Unit Hydrograph transfer function, defined
#' by Auto-Regressive and Moving Average coefficients.
#'
#' The transfer function used here, with input \var{u} and output \var{x} is:
#' \deqn{x[t] = a_1 x[t-1] + \ldots + a_n x[t-n] + }{ x[t] = a[1] x[t-1] + ...
#' + a[n] x[t-n] + b[0] u[t-d] + ... + b[m] u[t-m-d]}\deqn{ b_0 u[t-\delta] +
#' \ldots + b_m u[t-m-\delta]}{ x[t] = a[1] x[t-1] + ... + a[n] x[t-n] + b[0]
#' u[t-d] + ... + b[m] u[t-m-d]}
#'
#' and the \emph{order} is denoted \eqn{(n, m)}, with delay \eqn{\delta}{d}.
#'
#' @aliases armax armax.sim ssg.armax normalise.armax
#' @param U input time series.
#' @param a_1,a_2,a_3,b_0,b_1,b_2,b_3 ARMAX coefficients. Auto-regressive terms
#' begin with \code{a} and moving average terms begin with \code{b}. See
#' Details section.
#' @param pars the ARMAX coefficients as a named vector. If this is given, it
#' will over-ride the named parmameter arguments. Any number of terms can be
#' given here, it is not limited to the named arguments.
#' @param delay lag (dead time) between input and response, in time steps.
#' @param init initial values for the autoregressive filter.
#' @param na.action function to remove missing values, e.g.
#' \code{\link[=na.omit.ts]{na.omit}}.
#' @param epsilon values smaller than this will be set to zero.
#' @param return_components whether to return exponential component time
#' series.  If \code{TRUE}, the parameters will be converted to an exponential
#' components formulation, and passed to \code{\link{expuh.sim}}. This may fail
#' in some cases.
#' @param theta the parameters as a named vector.
#' @return the model output as a \code{\link{ts}} object, with the same
#' dimensions and time window as the input \code{U}.  If
#' \code{return_components = TRUE}, it will have multiple columns named
#' \code{Xs}, \code{Xq} and, if relevant, \code{X3}.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{armax.sriv.fit}}, \code{\link{arima}}
#' @references Jakeman, A.J., I.G. Littlewood, and P.G. Whitehead (1990),
#' Computation of the instantaneous unit hydrograph and identifiable component
#' flows with application to two small upland catchments, \emph{Journal of
#' Hydrology}, 117: 275-300.
#' @keywords ts
#' @examples
#'
#' data(HydroTestData)
#' fit <- hydromad(HydroTestData, routing = "armax",
#'                 rfit = list("ls", order = c(n = 2, m = 1)))
#' pars <- coef(fit)
#' pars
#'
#' xyplot(armax.sim(HydroTestData[,"P"], pars = pars))
#'


#' @importFrom stats filter na.pass
#'
#' @export
armax.sim <-
  function(U, a_1 = 0, a_2 = 0, a_3 = 0,
           b_0 = 1, b_1 = 0, b_2 = 0, b_3 = 0,
           pars = NULL,
           delay = 0, init = 0, na.action = na.pass,
           epsilon = hydromad.getOption("sim.epsilon"),
           return_components = FALSE) {
    ## parameter vectors
    a <- c(a_1, a_2, a_3)
    b <- c(b_0, b_1, b_2, b_3)
    a <- stripzeros(a)
    b <- stripzeros(b, up.to = 1)
    ## allow parameters to be passed in a single named vector
    if (length(pars) > 0) {
      pars <- tfParsConvert(pars, "a,b")
      a <- pars[grep("^a", names(pars))]
      b <- pars[grep("^b", names(pars))]
      stopifnot(length(b) > 0)
    }
    ## model order
    n <- length(a)
    m <- length(b) - 1
    ## check values
    stopifnot(all(is.finite(c(a, b))))
    ## stability check (based on stats::arima)
    if (length(a) > 0) {
      if (!all(Mod(polyroot(c(1, -a))) > .95)) {
        stop("AR component not in the region of stationarity")
      }
    }
    ## note U is allowed to be multi-variate, i.e. multiple columns
    U <- na.action(U)
    ## to return components, need to convert and pass to expuh.sim
    if (return_components) {
      abpars <- c(a, b)
      names(abpars) <- c(
        paste("a", seq_len(n), sep = "_"),
        paste("b", seq(0, m), sep = "_")
      )
      return(expuh.sim(U,
        pars = tfParsConvert(abpars, "tau,v"),
        delay = delay, Xs_0 = init[1],
        epsilon = epsilon,
        return_components = TRUE
      ))
    }
    inAttr <- attributes(U)
    U <- as.ts(U)
    ## apply 'U' delay in reverse to 'X' (i.e. lag by 'delay' steps)
    ## so can take delay as 0 for simulation purposes
    if (delay != 0) {
      U <- lag(U, -delay)
    }
    ## run ARMAX model
    init <- rep(init, length = n)
    X <- U
    X[] <- filter(U, b, sides = 1)
    if (length(a) > 0) {
      X <- shiftWindow(X, -m, fill = 0)
      X[] <- filter_ok(X, a, method = "recursive", init = init)
      X <- shiftWindow(X, m)
    }
    ## zap simulated values smaller than epsilon
    X[abs(X) < epsilon] <- 0
    X <- shiftWindow(X, delay)
    attributes(X) <- inAttr
    X
  }


#' @rdname armax.sim
#' @export
ssg.armax <- function(theta) {
  ssg.tf.coef(theta)
}

#' @rdname armax.sim
#' @export
normalise.armax <- function(theta) {
  normalise.tf.coef(theta)
}
