## Xs and X3 uses rolling mean effective rainfall as input
## only series=0 is allowed
maexpuh.sim <-
  function(U, delay = 0, tau_s = 0, tau_q = 0, tau_3 = 0, v_s = 1,
           v_q = NA, v_3 = 0, series = 0, loss = 0, Xs_0 = 0, Xq_0 = 0,
           X3_0 = 0, w_s = 7, w_3 = 30,
           pars = NULL, return_components = FALSE, na.action = na.pass,
           epsilon = hydromad.getOption("sim.epsilon")) {
    if (!is.null(pars)) {
      ccall <- match.call()
      ccall$pars <- NULL
      ccall <- as.call(modifyList(as.list(ccall), as.list(pars)))
      return(eval.parent(ccall))
    }
    delay <- round(delay)
    series <- round(series)
    stopifnot(series == 0)
    U <- na.action(U)
    if (delay != 0) {
      U <- lag(U, -delay)
    }
    if (is.na(v_q)) {
      if (series == 0) {
        v_q <- 1 - v_s - v_3
      }
      v_q <- max(0, min(1, v_q))
    }
    stopifnot(all(c(tau_s, tau_q, tau_3) >= 0))
    alpha_s <- exp(-1 / tau_s)
    alpha_q <- exp(-1 / tau_q)
    alpha_3 <- exp(-1 / tau_3)
    beta_s <- v_s * (1 - alpha_s)
    beta_q <- v_q * (1 - alpha_q)
    beta_3 <- v_3 * (1 - alpha_3)
    Xs <- Xq <- X3 <- Us <- U3 <- Xstemp <- X3temp <- U * NA
    lossVal <- (1 - alpha_s) * loss
    ## Pre and post filtered Xs
    Us[] <- filter(U, rep(1 / w_s, w_s), sides = 1)
    Us[1:(w_s - 1)] <- cumsum(U[1:(w_s - 1)]) / 1:(w_s - 1)
    Xstemp[] <- filter_loss(beta_s * U, alpha_s,
      loss = lossVal,
      init = Xs_0
    )
    Xs[] <- filter(Xstemp, rep(1 / w_s, w_s), sides = 1)
    Xs[1:(w_s - 1)] <- cumsum(Xstemp[1:(w_s - 1)]) / 1:(w_s - 1)
    if ((series == 0) || return_components) {
      Xq[] <- filter(beta_q * U, alpha_q,
        method = "recursive",
        init = Xq_0
      )
      if (v_3) {
        ## Pre and post filtered X3
        U3[] <- filter(U, rep(1 / w_3, w_3), sides = 1)
        U3[1:(w_3 - 1)] <- cumsum(U[1:(w_3 - 1)]) / 1:(w_3 - 1)
        X3temp[] <- filter(beta_3 * U3, alpha_3,
          method = "recursive",
          init = X3_0
        )
        X3[] <- filter(X3temp, rep(1 / w_3, w_3), sides = 1)
        X3[1:(w_3 - 1)] <- cumsum(X3temp[1:(w_3 - 1)]) / 1:(w_3 - 1)
      }
    }
    if (v_3 == 0) {
      X3 <- NULL
    }
    Xs <- shiftWindow(Xs, delay)
    Xq <- shiftWindow(Xq, delay)
    X3 <- if (!is.null(X3)) {
      shiftWindow(X3, delay)
    }
    Xs[Xs < epsilon] <- 0
    Xq[Xq < epsilon] <- 0
    if (!is.null(X3)) {
      X3[X3 < epsilon] <- 0
    }
    if (return_components) {
      if (v_3) {
        return(cbind(Xs = Xs, Xq = Xq, X3 = X3))
      } else {
        return(cbind(Xs = Xs, Xq = Xq))
      }
    }
    else {
      if (!is.null(X3)) {
        return(Xs + Xq + X3)
      } else {
        return(Xs + Xq)
      }
    }
  }
