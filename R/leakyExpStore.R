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
