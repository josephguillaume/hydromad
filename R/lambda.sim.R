## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

#' Transfer function with two exponential components and variable partitioning
#'
#' \var{Lambda} unit hydrograph.  Transfer function with two exponential
#' components and variable partitioning.
#'
#'
#' The \var{lambda} unit hydrograph model is a variant of the second-order
#' \code{\link{expuh}} model, i.e. two exponentially receding stores in
#' parallel. The \var{lambda} form allows the partitioning of flow between
#' quick and slow components to depend on the magnitude of effective rainfall.
#' In this model, runoff from large rainfall events tends to be quick flow, and
#' runoff from small events tends to be slow flow.
#'
#' \deqn{v_s[t] = v_{s,0} U[t] ^ \lambda} \deqn{v_q[t] = 1 - v_s[t]}
#'
#' where \var{U} is the input (effective rainfall); \eqn{v_{s,0}} is the
#' maximum fractional volume of the slow flow component, and is given by the
#' \code{v_s} argument.
#'
#' The \eqn{\lambda} parameter (\code{lambda} argument) must be between 0 and
#' -1; the case \code{lambda = 0} corresponds to the basic \code{\link{expuh}}
#' model.
#'
#' @name lambda
#' @aliases lambda.sim lambda.inverse.sim lambda.inverse.fit normalise.lambda ssg.lambda
#' @param U input time series.
#' @param delay lag (dead time) between input and response, in time steps.
#' @param tau_s,tau_q time constants for the exponential components.
#' @param lambda variable partitioning parameter, see Details.
#' @param v_s maximum fractional volume for the slower exponential component.
#' @param loss a constant loss (or gain) term subtracted from the \emph{slow}
#' (\code{s}) component.
#' @param Xs_0,Xq_0 initial values of the exponential components.
#' @param return_components whether to return all component time series.
#' @param na.action function to remove missing values, e.g.
#' \code{\link[=na.omit.ts]{na.omit}}.
#' @param epsilon values smaller than this will be set to zero.
#' @return the model output as a \code{\link{ts}} object, with the same
#' dimensions and time window as the input \code{U}.  If
#' \code{return_components = TRUE}, it will have multiple columns named
#' \code{Xs} and \code{Xq}.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{expuh}}, \code{\link{lambda.inverse.sim}}
#' @references ...
#' @keywords ts
#'
#' @export
lambda.sim <-
  function(U, delay = 0,
           tau_s = 0, tau_q = 0,
           lambda = 0, v_s = 1,
           loss = 0,
           Xs_0 = 0, Xq_0 = 0,
           return_components = FALSE,
           na.action = na.pass,
           epsilon = hydromad.getOption("sim.epsilon")) {
    delay <- round(delay)
    ## note U is allowed to be multi-variate, i.e. multiple columns
    U <- na.action(U)
    ## apply 'U' delay in reverse to 'X' (i.e. lag by 'delay' steps)
    ## so can take delay as 0 for simulation purposes
    if (delay != 0) {
      U <- lag(U, -delay)
    }
    ## check values
    stopifnot(all(c(tau_s, tau_q) >= 0))
    stopifnot((-1 <= lambda) && (lambda <= 0))
    ## convert from 'tau' and 'v' to 'alpha' and 'beta'
    alpha_s <- exp(-1 / tau_s)
    alpha_q <- exp(-1 / tau_q)
    ## lambda parameter defines dependence of v_s on U
    v_s <- v_s * (U^lambda)
    v_s <- pmax(pmin(v_s, 1), 0) ## ensure (0 <= v_s <= 1)
    v_q <- pmax(1 - v_s, 0)
    ## note: here v_s / v_q and beta_s / beta_q are vectors!
    beta_s <- v_s * (1 - alpha_s)
    beta_q <- v_q * (1 - alpha_q)
    ## apply exponential decay filter to each component
    ## note filter_loss is equivalent to filter if loss = 0
    ## convert loss from G[k] model into a Q[k] formulation
    lossVal <- (1 - alpha_s) * loss
    Xs <- filter_loss(beta_s * U, alpha_s, loss = lossVal, init = Xs_0)
    Xq <- beta_q * U
    Xq[] <- filter_ok(Xq, alpha_q, method = "recursive", init = Xq_0)

    ## align results to original input
    Xs <- shiftWindow(Xs, delay)
    Xq <- shiftWindow(Xq, delay)

    ## zap simulated values smaller than epsilon
    Xs[Xs < epsilon] <- 0
    Xq[Xq < epsilon] <- 0

    if (return_components) {
      return(cbind(Xs = Xs, Xq = Xq))
    } else {
      return(Xs + Xq)
    }
  }

#' @export
ssg.lambda <- function(theta) {
  return(1)
}

#' @export
normalise.lambda <- function(theta) {
  return(theta)
}
