## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

expuh.sim <-
  function(U, delay = 0,
           tau_s = 0, tau_q = 0, tau_3 = 0,
           v_s = 1, v_q = NA, v_3 = 0,
           series = 0, loss = 0,
           Xs_0 = 0, Xq_0 = 0, X3_0 = 0,
           pars = NULL,
           return_components = FALSE,
           na.action = na.pass,
           epsilon = hydromad.getOption("sim.epsilon")) {
    if (!is.null(pars)) {
      ccall <- match.call()
      ccall$pars <- NULL
      ccall <- as.call(modifyList(
        as.list(ccall),
        as.list(pars)
      ))
      return(eval.parent(ccall))
    }
    delay <- round(delay)
    series <- round(series)
    stopifnot(series %in% 0:3)
    if ((series > 1) && (v_3 == 0)) {
      stop("series > 1 is only valid with v_3 > 0")
    }
    ## note U is allowed to be multi-variate, i.e. multiple columns
    U <- na.action(U)
    ## apply 'U' delay in reverse to 'X' (i.e. lag by 'delay' steps)
    ## so can take delay as 0 for simulation purposes
    if (delay != 0) {
      U <- lag(U, -delay)
    }

    ## default value of v_q depends on 'series'
    ## (by convention only, this defines configuration of the system)
    if (is.na(v_q)) {
      if (series == 0) {
        ## all in parallel
        v_q <- 1 - v_s - v_3
      }
      if (series == 1) {
        v_q <- 1
      }
      if (series == 2) {
        v_q <- 1 - v_s
      }
      if (series == 3) {
        v_q <- 1
      }
      v_q <- max(0, min(1, v_q))
    }

    stopifnot(all(c(tau_s, tau_q, tau_3) >= 0))
    ## convert from 'tau' and 'v' to 'alpha' and 'beta'
    alpha_s <- exp(-1 / tau_s)
    alpha_q <- exp(-1 / tau_q)
    alpha_3 <- exp(-1 / tau_3)
    beta_s <- v_s * (1 - alpha_s)
    beta_q <- v_q * (1 - alpha_q)
    beta_3 <- v_3 * (1 - alpha_3)
    ## apply exponential decay filter to each flow component.
    ## components will be added together at the end.
    Xs <- Xq <- X3 <- U * NA
    ## default ('slow') flow component optionally includes loss term.
    ## note filter_loss is equivalent to filter if loss = 0
    ## convert loss from G[k] model into a Q[k] formulation
    lossVal <- (1 - alpha_s) * loss
    Xs[] <- filter_loss(beta_s * U, alpha_s, loss = lossVal, init = Xs_0)

    if ((series == 0) || return_components) {
      ## components in parallel.
      Xq[] <- filter_ok(beta_q * U, alpha_q, method = "recursive", init = Xq_0)
      if (v_3) {
        X3[] <- filter_ok(beta_3 * U, alpha_3, method = "recursive", init = X3_0)
      }
    } else {
      if (v_3 == 0) {
        ## second-order model in series (s --> q)
        Xq[] <- filter_ok(beta_q * Xs, alpha_q, method = "recursive", init = Xq_0)
        Xs[] <- 0
      } else {
        ## third-order model
        if (series == 1) {
          ## two components in series and one in parallel
          ## (q & 3 are in series; s in parallel)
          X3[] <- filter_ok(beta_3 * U, alpha_3, method = "recursive", init = X3_0)
          Xq[] <- filter_ok(beta_q * X3, alpha_q, method = "recursive", init = Xq_0)
          X3[] <- 0
        } else if (series == 2) {
          ## one component in series with two in parallel
          ## (3 in series; s & q in parallel)
          Xq[] <- filter_ok(beta_q * U, alpha_q, method = "recursive", init = Xq_0)
          Xq[] <- filter_ok(beta_3 * Xq, alpha_3, method = "recursive", init = X3_0)
          Xs[] <- filter_ok(beta_3 * Xs, alpha_3, method = "recursive", init = X3_0)
          X3[] <- 0
        } else if (series == 3) {
          ## three components in series
          Xq[] <- filter_ok(beta_q * Xs, alpha_q, method = "recursive", init = Xq_0)
          X3[] <- filter_ok(beta_3 * Xq, alpha_3, method = "recursive", init = X3_0)
          Xs[] <- 0
          Xq[] <- 0
        } else {
          stop("unrecognised values of 'series': ", series)
        }
      }
    }

    if (v_3 == 0) X3 <- NULL

    ## align results to original input
    Xs <- shiftWindow(Xs, delay)
    Xq <- shiftWindow(Xq, delay)
    X3 <- if (!is.null(X3)) shiftWindow(X3, delay)

    ## zap simulated values smaller than epsilon
    Xs[Xs < epsilon] <- 0
    Xq[Xq < epsilon] <- 0
    if (!is.null(X3)) X3[X3 < epsilon] <- 0

    if (return_components) {
      if (v_3) {
        return(cbind(Xs = Xs, Xq = Xq, X3 = X3))
      } else {
        return(cbind(Xs = Xs, Xq = Xq))
      }
    } else {
      if (!is.null(X3)) {
        return(Xs + Xq + X3)
      } else {
        return(Xs + Xq)
      }
    }
  }


ssg.expuh <- function(theta) {
  if (length(theta) == 0) {
    return(1)
  }
  theta <- tfParsConvert(theta, "a,b")
  ssg.tf.coef(theta)
}

normalise.expuh <- function(theta) {
  theta <- tfParsConvert(theta, "a,b")
  tmp <- normalise.tf.coef(theta)
  tfParsConvert(tmp, "tau,v")
}
