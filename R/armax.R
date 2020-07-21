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
