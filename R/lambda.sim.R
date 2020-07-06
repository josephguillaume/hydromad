## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

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


ssg.lambda <- function(theta) {
  return(1)
}

normalise.lambda <- function(theta) {
  return(theta)
}
