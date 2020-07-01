## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


## TODO: keep attributes of Q

lambda.inverse.sim <-
  function(Q, P = NULL,
           pars = c(tau_s = 0, tau_q = 0, v_s = 1, v_3 = 0),
           delay = 0, Xs_0 = 0, Xq_0 = 0,
           rain.factor = 1.1,
           rises.only = FALSE,
           use.Qm = TRUE,
           mass.balance = FALSE,
           scale.window = NA) {
    if (!is.ts(Q)) Q <- as.ts(Q)
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
        P <- shiftWindow(P, delay = delay)
        ## (and fix up time series attributes)
        P <- lag(P, -delay)
      }
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


    ###### TODO
    if (TRUE || hydromad.getOption("pure.R.code")) {
      ## slow version in R for cross-checking
      alpha_s <- exp(-1 / pars[["tau_s"]])
      alpha_q <- exp(-1 / pars[["tau_q"]])
      if (is.na(alpha_q)) alpha_q <- 0
      v_s_0 <- pars[["v_s"]]
      lambda <- pars[["lambda"]]

      Qs_0 <- min(Q[1:10])
      Qq_0 <- 0
      Qst1 <- Qs_0
      Qqt1 <- Qq_0
      cDiffFun <- function(U, Qstar) {
        Qtest <- (1 - alpha_q) * U +
          v_s_0 * (alpha_q - alpha_s) * U^(lambda + 1)
        return(abs(Qstar - Qtest))
      }
      for (t in 1:length(Q)) {
        Qstar <- Q[t] - alpha_s * Qst1 - alpha_q * Qqt1
        ## trivial case: P is zero, therefore U is zero
        if ((P[t] == 0) || (Qstar <= 0)) {
          U[t] <- 0
        } else {
          ## estimate U[t] from Q[t], Q[t-1], ...
          ## we need to know v_s, but this depends on U[t] ^ lambda.
          ## constraints:
          v_s_max <- min(1, v_s_0)
          ## bound based on v_s as function of U (U <= P)
          v_s_min <- max(0, min(v_s_max, v_s_0 * P[t]^lambda))
          # v_s_max <- max(v_s_max, v_s_min)
          beta_max <- v_s_min * (1 - alpha_s) + (1 - v_s_min) * (1 - alpha_q)
          beta_min <- v_s_max * (1 - alpha_s) + (1 - v_s_max) * (1 - alpha_q)
          Umin <- Qstar / beta_max
          Umax <- Qstar / beta_min
          ## bound based on rainfall
          Umax <- min(Umax, P[t])
          ## bound based on lambda = -1
          Umin <- max(Umin, (Qstar - v_s_0 * (alpha_q - alpha_s)) /
            (1 - alpha_q))
          Umin <- min(Umin, Umax)
          ## ensure that U values imply possible values of v_s
          if ((v_s_0 * Umax^lambda) >= 1) {
            Umin <- Umax
          }
          ## check whether we have constrained U to a unique solution
          if ((Umax - Umin) < 0.0001) {
            U[t] <- Umax
          } else {
            ## approximation based on linear interpolation - any use?
            # Uest <- Umin + (Umax - Umin) * (1+lambda)
            result <- optimise(cDiffFun, interval = c(Umin, Umax), Qstar = Qstar)
            U[t] <- result$minimum
          }
        }
        ## calculate flow components
        v_s_t <- 0
        if (U[t] > 0) {
          v_s_t <- v_s_0 * (U[t]^lambda)
        }
        v_s_t <- max(0, min(1, v_s_t))
        beta_s_t <- v_s_t * (1 - alpha_s)
        beta_q_t <- (1 - v_s_t) * (1 - alpha_q)
        Qst1 <- alpha_s * Qst1 + beta_s_t * U[t]
        Qqt1 <- alpha_q * Qqt1 + beta_q_t * U[t]
      }
    } else {
      Qs_0 <- min(Q[1:10])
      Qq_0 <- 0
      U <- .C(inverse_filter_lambda,
        as.double(Q),
        as.integer(length(Q)),
        as.double(Qs_0),
        as.double(Qq_0),
        as.double(pars["tau_s"]),
        as.double(pars["tau_q"]),
        as.double(pars["v_s"]),
        as.double(pars["lambda"]),
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
      U <- shiftWindow(U, delay = -delay)
      ## (and fix up time series attributes)
      U <- lag(U, delay)
    }
    U
  }
