## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Joseph Guillaume <joseph.guillaume@anu.edu.au>
##

expuh3s.sim <-
  function(U, delay = 0, v_s,
           tau_q, tau_s, tau_g,
           R, G_1, loss, G_2,
           Xs_0 = 0, Xq_0 = 0, Xg_0 = 0,
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
    ## note U is allowed to be multi-variate, i.e. multiple columns
    U <- na.action(U)
    ## apply 'U' delay in reverse to 'X' (i.e. lag by 'delay' steps)
    ## so can take delay as 0 for simulation purposes
    if (delay != 0) {
      U <- lag(U, -delay)
    }

    stopifnot(all(c(tau_s, tau_q, tau_g) >= 0))

    ## convert from 'tau' and 'v' to 'alpha' and 'beta'
    alpha_q <- exp(-1 / tau_q)
    beta_q <- (1 - v_s) * (1 - alpha_q)

    ### apply exponential decay filter to quick flow
    Xq <- U * NA
    Xq[] <- filter_ok(beta_q * U, alpha_q, method = "recursive", init = Xq_0)

    ## Upper store
    upper <- leakyExpStore.sim(v_s * U, tau_s, loss = R, thres = G_1, init = Xs_0, return_components = TRUE)
    Xs <- upper[, "Q"]

    ## Lower store
    lower <- leakyExpStore.sim(upper[, "L"], tau_g,
      loss = loss, thres = G_2, init = Xg_0,
      return_components = FALSE
    )
    Xg <- lower[, "Q"]

    ## align results to original input
    Xs <- shiftWindow(Xs, delay)
    Xq <- shiftWindow(Xq, delay)
    Xg <- shiftWindow(Xg, delay)

    ## zap simulated values smaller than epsilon
    Xs[Xs < epsilon] <- 0
    Xq[Xq < epsilon] <- 0
    Xg[Xg < epsilon] <- 0

    ## TODO: return store as well as flows
    if (return_components) {
      return(cbind(Xs = Xs, Xq = Xq, Xg = Xg))
    } else {
      return(Xs + Xq + Xg)
    }
  }
