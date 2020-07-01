## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

expuh.inverse.sim <-
  function(DATA, delay = 0,
           tau_s = 0, tau_q = 0, tau_3 = 0,
           v_s = 1, v_q = NA, v_3 = 0,
           series = 0,
           Xs_0 = 0, Xq_0 = 0, X3_0 = 0,
           pars = NULL,
           ...) {
    pars0 <- list(
      tau_s = tau_s, tau_q = tau_q, tau_3 = tau_3,
      v_s = v_s, v_q = v_q, v_3 = v_3, series = series
    )
    pars <- modifyList(pars0, as.list(pars))
    pars <- tfParsConvert(pars, "a,b")
    armax.inverse.sim(DATA, pars = pars, ...)
  }


armax.inverse.sim <-
  function(DATA,
           a_1 = 0, a_2 = 0, a_3 = 0,
           b_0 = 1, b_1 = 0, b_2 = 0, b_3 = 0,
           pars = NULL,
           delay = 0, init = 0,
           rain.factor = 1.1,
           rises.only = FALSE,
           use.Qm = TRUE,
           use.fft.method = FALSE,
           constrain.fft = TRUE,
           mass.balance = use.fft.method,
           scale.window = NA) {
    P <- NULL
    if (NCOL(DATA) > 1) {
      Q <- DATA[, "Q"]
      if ("P" %in% colnames(DATA)) {
        P <- DATA[, "P"]
      }
    } else {
      Q <- DATA
    }
    inAttr <- attributes(Q)
    Q <- as.ts(Q)
    if (!is.null(P)) P <- as.ts(P)

    pars <- tfParsConvert(pars, "a,b")
    pars0 <- list(
      a_1 = a_1, a_2 = a_2, a_3 = a_3,
      b_0 = b_0, b_1 = b_1, b_2 = b_2, b_3 = b_3
    )
    pars <- unlist(modifyList(pars0, as.list(pars)))
    tfParsCheck(pars)
    ## extract first parts of parameter names (before underscore)
    parSymbols <- gsub("_.*", "", names(pars))
    a <- pars[parSymbols == "a"]
    b <- pars[parSymbols == "b"]
    a <- stripzeros(a)
    b <- stripzeros(b, up.to = 1)
    n <- length(a)
    m <- length(b) - 1
    # if (constrain.zeros) {
    #    if (!is.null(P)) {
    #        zeros <- mean(zapsmall(P, 2) == 0, na.rm = TRUE)
    #    } else {
    #        zeros <- mean(zapsmall(diff(Q), 2) <= 0, na.rm = TRUE)
    #    }
    # }
    ## set NAs in flow to zero (insert them again afterwards)
    isna <- is.na(Q)
    Q[isna] <- 0
    ## initialise P
    if (!is.null(P)) {
      if (!is.ts(P)) P <- as.ts(P)
      P <- P * rain.factor
      ## apply delay
      if (delay > 0) {
        P <- shiftWindow(P, delay = delay, and.lag = TRUE)
      }
      P <- window(P, start = start(Q), end = end(Q), extend = TRUE)
    } else {
      P <- Q
      P[] <- 1e12 # Inf
    }
    if (rises.only) {
      ## constrain U to time steps with rising flow
      diffrising <- (diff(Q) > 0)
      ## backwards difference, i.e. rises are placed at their end time
      rising <- c(FALSE, diffrising)
      P[!rising] <- 0
    }
    ## restrict flow on missing rain days
    P[is.na(P)] <- 0
    stopifnot(length(P) == length(Q))
    ## initialise U
    U <- Q
    U[] <- 0

    if (use.fft.method) {
      ## fourier deconvolution method
      bb <- b[-1] / b[1]
      maform <- ARMAtoMA(a, bb, lag.max = length(Q))
      U <- wiener(Q, maform, N = 0)
      U <- pmax(U, 0)
      # U[U <= quantile(U, zeros, na.rm = TRUE)] <- 0
      if (constrain.fft) {
        U <- pmin(U, P)
      }
    } else

    if (hydromad.getOption("pure.R.code")) {
      ## slow version in R for cross-checking
      Qm <- Q
      for (t in seq(max(n, m) + 1, length(Q))) {
        Ut <- Q[t]
        if (n > 0) {
          Ut <- Ut - sum(a * Qm[t - seq_len(n)])
        }
        if (m > 0) {
          Ut <- Ut - sum(b[-1] * U[t - seq_len(m)])
        }
        Ut <- Ut / b[1]
        U[t] <- max(0, min(Ut, P[t]))
        if (use.Qm) {
          Qm[t] <- (sum(a * Qm[t - seq_len(n)]) +
            sum(b * U[t - 0:m]))
        }
      }
    } else {
      U <- .C(inverse_filter,
        as.double(Q[]), ## force copy
        as.integer(length(Q)),
        as.double(a),
        as.double(b),
        as.integer(n),
        as.integer(m),
        as.integer(use.Qm),
        as.double(P),
        U = double(length(Q)),
        DUP = FALSE, PACKAGE = "hydromad"
      )$U
      ## make it a time series object again
      attributes(U) <- attributes(Q)
    }
    U[isna] <- NA
    ## scale to ensure mass balance with Q
    if (mass.balance) {
      if (is.na(scale.window)) {
        U <- U * sum(Q[!isna]) / sum(U[!isna])
      } else {
        sc <- simpleSmoothTs(cbind(Q, U), width = scale.window)
        sc <- sc[, "Q"] / sc[, "U"]
        sc <- na.locf(na.locf(sc, na.rm = FALSE), fromLast = TRUE, na.rm = FALSE)
        U <- U * sc
      }
    }
    ## inverse.sim should be in sync with P not Q,
    ## so delay still applies
    if (delay > 0) {
      U <- shiftWindow(U, delay = -delay, and.lag = TRUE)
    }
    attributes(U) <- inAttr
    U
  }


## based on http://en.wikipedia.org/w/index.php?title=Wiener_filter&oldid=211541130
wiener <- function(y, h, N, gamma = 0) {
  length(h) <- length(y)
  h[is.na(h)] <- 0
  y[is.na(y)] <- 0
  H <- fft(h)
  Y <- fft(y)
  S <- abs(Y)^2 / length(y)^2

  ## experimental threshold to reduce sensitivity
  ## based on http://cnx.org/content/m13144/latest/wienerFilter.m
  ## TODO: is this ever useful?
  zap <- Mod(H) < gamma
  if (any(zap)) {
    H[zap] <- gamma * Mod(H[zap]) / H[zap]
  }

  G <- (Conj(H) * S) /
    (abs(H)^2 * S + N)

  x <- Re(fft(G * Y, inverse = TRUE) / length(y))
  x
}
