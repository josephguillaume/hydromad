#' Exponential store with zero flows and thresholded loss
#'
#' An exponential store (linear transfer function) which has a loss term,
#' produces no flow when the store drops below a level, and can therefore model
#' longer-term disconnection of a store from streamflow.
#'
#' Storage is increased by effective rainfall and decreased by flow and losses.
#' \deqn{G[k] = G[k-1] + U[k] - Q[k] - L[k]}
#'
#' Flow is proportional to storage \deqn{Q[k] = a G[k], G[k]>0} \deqn{Q[k] = 0,
#' otherwise}
#'
#' Loss switches off at some threshold, as a piece-wise continuous function.
#' \deqn{L[k] = L, G[k]>T+L (with T<=0)} \deqn{L[k] = G[k] - T, T+L > G[k] > T}
#' \deqn{L[k] = 0, G[k]<T}
#'
#' @name leakyExpStore
#' @aliases leakyExpStore leakyExpStore.sim
#' @param x input time series.
#' @param tau Time constant for flow from exponential store.
#' @param loss Constant loss that occurs while the value of the store is
#' greater than \code{thres}. See details.
#' @param thres Trigger level to turn off loss, which can be interpreted as a
#' store capacity.
#' @param init Initial value of exponential store
#' @param return_components whether to return store value (G) as well as flow
#' (Q)
#' @return the flow time series with the same dimensions and time windows as
#' the input \code{x}.  If \code{return_components = TRUE}, it will have
#' multiple columns named \code{G} and \code{Q}.
#' @author Joseph Guillaume \email{joseph.guillaume@@anu.edu.au} with advice
#' from Barry Croke
#' @seealso \code{\link{expuh3s.sim}}
#' @keywords ts
#' @examples
#'
#' U <- ts(c(1, rep(0, 10), 1, rep(0, 20)))
#'
#' ## Without a loss, equivalent to expuh
#' ##  Loss threshold has no effect
#' all.equal(leakyExpStore.sim(U, 5, loss = 0, thres = 0), expuh.sim(U, tau_s = 5))
#'
#' ## With losses stopping when flow stops, equivalent to expuh with a loss
#' all.equal(
#'   leakyExpStore.sim(U, 5, loss = 0.1, thres = 0),
#'   expuh.sim(U, tau_s = 5, loss = 0.1)
#' )
#'
#' ## Plot of unit hydrographs
#' xyplot(cbind(
#'   "U" = U,
#'   "loss=0" = expuh.sim(U, tau_s = 5),
#'   "loss=0.1 & thres=0" = expuh.sim(U, tau_s = 5, loss = 0.1),
#'   "loss=0.1 & thres=-0.3" = leakyExpStore.sim(U, 5, loss = 0.1, thres = -0.3),
#'   "loss=0.1 & thres=-Inf" = leakyExpStore.sim(U, 5, loss = 0.1, thres = -Inf)
#' ), superpose = TRUE)
#'
#' ## Time series plot of value of store
#' xyplot(cbind(
#'   "thres=0" = leakyExpStore.sim(U, 5,
#'     loss = 0.1, thres = 0,
#'     return_components = TRUE
#'   )[, "G"],
#'   "thres=-0.3" = leakyExpStore.sim(U, 5,
#'     loss = 0.1, thres = -0.3,
#'     return_components = TRUE
#'   )[, "G"],
#'   "thres=-Inf" = leakyExpStore.sim(U, 5,
#'     loss = 0.1, thres = -Inf,
#'     return_components = TRUE
#'   )[, "G"]
#' ), superpose = TRUE)
#'
#' ## Time series of loss
#' xyplot(cbind(
#'   "thres=0" = leakyExpStore.sim(U, 5,
#'     loss = 0.1, thres = 0,
#'     return_components = TRUE
#'   )[, "L"],
#'   "thres=-0.3" = leakyExpStore.sim(U, 5,
#'     loss = 0.1, thres = -0.3,
#'     return_components = TRUE
#'   )[, "L"],
#'   "thres=-Inf" = leakyExpStore.sim(U, 5,
#'     loss = 0.1, thres = -Inf,
#'     return_components = TRUE
#'   )[, "L"]
#' ), superpose = TRUE)
#' @export
leakyExpStore.sim <- function(x, tau, loss, thres, init = 0, return_components = FALSE) {
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(loss))
  stopifnot(is.numeric(thres))
  stopifnot(is.numeric(init))
  stopifnot(is.numeric(tau))
  stopifnot(length(tau) == 1)
  a <- exp(1 / tau) - 1
  ## stopifnot(length(x) > length(a) + 1) FIXME
  xAttrs <- attributes(x)
  bad <- !is.finite(x)
  x[bad] <- 0
  stopifnot(thres <= 0)
  G <- c(init, x * 0)
  Q <- c(NA, x * 0)
  x <- c(init, x)
  L <- rep(loss, length(x))

  ## COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
  COMPILED <- FALSE
  if (COMPILED) {
    ## TODO
    ## ans <- .C(leakyExpStore
  } else {
    for (k in 2:length(x)) {
      G0 <- G[k - 1] + x[k] - loss
      if (G0 <= 0) {
        G[k] <- G0
        Q[k] <- 0
      } else if (G0 > 0) {
        G[k] <- G0 / (1 + a)
        Q[k] <- a * G[k]
      }
      ## loss only switches off after/when Q switches off
      if (G[k] < thres) {
        L[k] <- L[k] + G[k] - thres
        G[k] <- thres
      }
    }
  }
  Q <- Q[-1]
  attributes(Q) <- xAttrs
  if (return_components) {
    G <- G[-1]
    L <- L[-1]
    attributes(G) <- xAttrs
    attributes(L) <- xAttrs
    return(cbind(G = G, Q = Q, L = L))
  } else {
    return(Q)
  }
}
