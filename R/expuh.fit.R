## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


#' Exponential components transfer function models
#'
#' A unit hydrograph (linear transfer function) defined as a system of
#' exponentially receding components. Each component is defined by its time
#' constant and fractional volume, and if there are multiple (up to 3) such
#' components they may be in a parallel and/or series configuration.
#'
#'
#' The \code{expuh} model is a transfer function translating an input time
#' series \var{U} into an output series \var{X}.  It describes a configuration
#' of exponentially decaying components, each defined by a recession rate
#' \eqn{\alpha} and peak response \eqn{\beta}. However, in hydrology these
#' parameters are more easily interpreted in terms of time constants \eqn{\tau}
#' (number of time steps to reduce to a fraction \eqn{1/e}, 37\%) and
#' fractional volumes \var{v}. These are directly related as:
#'
#' \deqn{\tau = -1 / \log(\alpha)}
#'
#' \deqn{v = \beta / (1 - \alpha)}
#'
#' If there are two components in parallel, these are conventionally called
#' slow (\var{s}) and quick (\var{q}) flow components. The total simulated flow
#' \var{X} is the sum of these; \eqn{X[t] = X_s[t] + X_q[t]}, and:
#'
#' \deqn{X_s[t] = \alpha_s X_s[t-1] + \beta_s U[t]} \deqn{X_q[t] = \alpha_q
#' X_q[t-1] + \beta_q U[t]}
#'
#' Two components might also be arranged in series rather than parallel, in
#' which case:
#'
#' \deqn{X_s[t] = \alpha_s X_s[t-1] + \beta_s U[t]} \deqn{X[t] = \alpha_q
#' X[t-1] + \beta_q X_s[t]}
#'
#' This configuration is specified by the argument \code{series = 1}. The
#' default \code{series = 0} specifies all components to be in parallel.
#'
#' In the case of three components, with corresponding time constants
#' \eqn{\tau_s}, \eqn{\tau_q} and \eqn{tau_3} (\code{tau_s, tau_q, tau_3}),
#' there are four possible types of configuration:
#'
#' \describe{ \item{list("series = 0")}{ all 3 components in parallel, i.e.
#' independent flows: \var{X = s + q + 3}. In this case \code{v_q} defaults to
#' \code{1 - v_s - v_3} in order to ensure that the total volume is 1.  }
#' \item{list("series = 1")}{ one component in parallel with two in series: the
#' \code{q} component is in series with the \code{3} component, and the
#' \code{s} component is in parallel: \var{X = (q * 3) + s}. In this case
#' \code{v_q} defaults to 1.  } \item{list("series = 2")}{ two components in
#' parallel with one in series: the \code{s} and \code{q} components are in
#' parallel and the \code{3} component is in series: \var{X = 3 * (s + q)}. In
#' this case \code{v_q} defaults to \code{1 - v_s} in order to ensure that the
#' total volume of the parallel component is 1. The total volume will be 1 if
#' \code{v_3} is also 1.  } \item{list("series = 3")}{ all 3 components in
#' series: \var{X = s * q * 3}.  In this case \code{v_q} defaults to 1. The
#' total volume will be 1 if \code{v_s} and \code{v_3} are also 1.  } }
#'
#' @importFrom stats optim
#'
#' @name expuh
#' @aliases ssg.expuh normalise.expuh expuh.ls.fit expuh.sim
#' @param DATA Placeholder
#' @param order Placeholder
#' @param quiet Placeholder
#' @param delay Placeholder
#' @param ... Placeholder
#' @return the model output as a \code{\link{ts}} object, with the same
#' dimensions and time window as the input \code{U}.  If
#' \code{return_components = TRUE}, it will have multiple columns named
#' \code{Xs}, \code{Xq} and, if relevant, \code{X3}.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{expuh.sriv.fit}} \code{\link{armax}}
#' @references Jakeman, A.J., I.G. Littlewood, and P.G. Whitehead (1990),
#' Computation of the instantaneous unit hydrograph and identifiable component
#' flows with application to two small upland catchments, \emph{Journal of
#' Hydrology}, 117: 275-300.
#' @keywords ts
#' @examples
#'
#' data(HydroTestData)
#' mod1 <- hydromad(HydroTestData,
#'   routing = "expuh",
#'   tau_s = 30, tau_q = 5, v_s = 0.5
#' )
#' flowcomps <- predict(mod1, return_components = TRUE)
#' xyplot(cbind(
#'   `Slow component` = flowcomps[, "Xs"],
#'   `Total flow` = flowcomps[, 1] + flowcomps[, 2]
#' ),
#' superpose = TRUE
#' ) +
#'   latticeExtra::layer(panel.refline(h = 0))
#'
#' U <- ts(c(1, rep(0, 30)))
#' xyplot(cbind(
#'   "tau_s = 10" = expuh.sim(U, tau_s = 10),
#'   "& tau_q = 1" = expuh.sim(U, tau_s = 10, tau_q = 1, v_s = 0.5),
#'   "&& v_s = 0.9" = expuh.sim(U, tau_s = 10, tau_q = 1, v_s = 0.9)
#' ),
#' superpose = TRUE
#' )
#' @export
expuh.ls.fit <-
  function(DATA,
           order = hydromad.getOption("order"),
           delay = hydromad.getOption("delay"),
           quiet = FALSE,
           ...) {
    model <- armax.ls.fit(DATA, order = order, delay = delay, ...)
    if (!inherits(model, "tf")) {
      return(model)
    }
    n <- order[1]
    poles <- arToPoles(coef(model)[seq_len(n)])
    eps <- sqrt(.Machine$double.eps)
    badpoles <- (Re(poles) < -eps) | (abs(Im(poles)) > eps)
    if (any(badpoles)) {
      if (!quiet) {
        message("armax fitted poles are non-physical; re-fitting with constraints")
      }
      model <-
        fitWithPoleConstraints(DATA,
          fitfun = armax.ls.fit, poles = poles,
          order = order, delay = delay, ...
        )
    }
    if (!inherits(model, "tf")) {
      return(model)
    }
    model$coefficients <- coef(model, "tau,v")
    model
  }


#' @export
expuh.sriv.fit <-
  function(DATA,
           order = hydromad.getOption("order"),
           delay = hydromad.getOption("delay"),
           quiet = FALSE,
           ...) {
    model <- armax.sriv.fit(DATA, order = order, delay = delay, ...)
    if (!inherits(model, "tf")) {
      return(model)
    }
    n <- order[1]
    poles <- arToPoles(coef(model)[seq_len(n)])
    eps <- sqrt(.Machine$double.eps)
    badpoles <- (Re(poles) < -eps) | (abs(Im(poles)) > eps)
    if (any(badpoles)) {
      if (!quiet) {
        message("armax fitted poles are non-physical; re-fitting with constraints")
      }
      model <-
        fitWithPoleConstraints(DATA,
          fitfun = armax.sriv.fit, poles = poles,
          order = order, delay = delay, ...
        )
    }
    if (!inherits(model, "tf")) {
      return(model)
    }
    model$coefficients <- coef(model, "tau,v")
    model
  }


#' @export
expuh.inverse.fit <-
  function(DATA,
           order = hydromad.getOption("order"),
           delay = hydromad.getOption("delay"),
           ...) {
    ## TODO: can do this directly?
    model <- armax.inverse.fit(DATA, order = order, delay = delay, ...)
    if (!inherits(model, "tf")) {
      return(model)
    }
    # model$vcov <- vcov(model)
    model$coefficients <- coef(model, "tau,v")
    model
  }



fitWithPoleConstraints <-
  function(DATA, fitfun, poles, ...,
           control = as.list(hydromad.getOption("optim.control.expuh"))) {
    ## non-physical roots; constrain them and fit again
    poles <- abs(Re(poles))
    poles <- pmax(poles, 1e-5)
    logpol0 <- log(poles)
    model <- structure("failed to fit expuh routing with constrained roots",
      class = "try-error"
    )
    bestVal <- Inf
    optFun <- function(logpol, ...) {
      pol <- exp(logpol)
      if (isTRUE(hydromad.getOption("catch.errors"))) {
        thisMod <- try(fitfun(DATA, ..., fixed.ar = polesToAr(pol)))
      } else {
        thisMod <- fitfun(DATA, ..., fixed.ar = polesToAr(pol))
      }
      if (!isValidModel(thisMod)) {
        return(NA)
      }
      val <- sum(abs(residuals(thisMod)), na.rm = TRUE)
      if (val < bestVal) {
        model <<- thisMod
      }
      val
    }
    if (!isTRUE(hydromad.getOption("catch.errors.optim"))) {
      try <- force
    } ## i.e. skip the try()
    # ans <- try(nlminb(logpol0, optFun, control = control, ...))
    ans <- try(optim(logpol0, optFun, control = control, ...))
    if (inherits(ans, "try-error")) {
      return(ans)
    }
    if (ans$convergence != 0) {
      msg <- paste("While re-fitting non-physical poles:", toString(ans$msg))
      if (!isTRUE(hydromad.getOption("quiet"))) {
        warning(msg)
      }
      model$msg <- msg
    }
    model
  }
