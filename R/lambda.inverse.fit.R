## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


lambda.inverse.fit <-
  function(DATA,
           order = c(n = 2, m = 1),
           delay = hydromad.getOption("delay"),
           fit.method = hydromad.getOption("inverse.fit.method"),
           normalise = NA,
           init.U = TRUE,
           pars = NULL,
           use.Qm = TRUE,
           rises.only = FALSE,
           ...,
           max.iterations = hydromad.getOption("inverse.iterations"),
           rel.tolerance = hydromad.getOption("inverse.rel.tolerance"),
           par.epsilon = hydromad.getOption("inverse.par.epsilon"),
           init.attempt = 0,
           trace = hydromad.getOption("trace")) {
    DATA <- as.ts(DATA)
    if ("Q" %in% colnames(DATA)) {
      Q <- DATA[, "Q"]
    } else if (NCOL(DATA) == 1) {
      Q <- DATA
    } else {
      stop("need Q in DATA")
    }
    P <- if ("P" %in% colnames(DATA)) DATA[, "P"]

    if (is.na(delay)) {
      delay <- 0
      if ("P" %in% colnames(DATA)) {
        delay <- estimateDelay(DATA, plot = FALSE)
      }
    }
    if (!is.null(pars)) {
      pars <- tfParsConvert(pars, "a,b")
      tfParsCheck(pars)
    }
    ## set NAs in flow to zero (insert them again afterwards)
    isna <- is.na(Q)
    Q[isna] <- 0
    sumQ <- sum(Q)

    U <- NULL
    ## option to initialise from U (as rises(Q)), rather than pars.init
    if (init.U) {

      #        if (!is.null(P)) {
      #            zeros <- mean(zapsmall(P, 2) == 0, na.rm = TRUE)
      #        } else {
      #            zeros <- mean(zapsmall(diff(Q), 2) <= 0, na.rm = TRUE)
      #        }
      ## backwards difference, i.e. rises are placed at their end time
      rises <- function(x) {
        x[] <- c(0, pmax(diff(x), 0))
        x
      }
      ## estimate U as the rises in Q, scaled
      #        U <- rises(Q)
      #        U[U <= quantile(U, zeros, na.rm=TRUE)] <- 0
      ## scale to ensure mass balance with Q
      #        U <- U * sum(Q[!isna]) / sum(U[!isna])
      if (!is.null(P)) {
        U <- P
      } else {
        U <- rises(Q)
        ## positive innovations from an ar1 model:
        # ar1 <- coef(arima(Q, order = c(1,0,0), include.mean = FALSE))
        # U <- pmax(Q - ar1 * lag(Q), 0)
      }

      ## estimate U as scaled P, scaled in a moving window
      ## (runoff coefficient)
      scale.window <- max(autocorrTime(Q) * 1.5, 16)
      sc <- simpleSmoothTs(cbind(Q, U), width = scale.window, c = 2)
      sc <- sc[, "Q"] / sc[, "U"]
      sc <- na.locf(na.locf(sc, na.rm = FALSE), fromLast = TRUE, na.rm = FALSE)
      sc[!is.finite(sc)] <- 0
      U <- U * sc
    } else {
      ## generate starting parameters
      if (is.null(pars)) {
        pars <- tf.pars.init(DATA,
          order = order, delay = delay,
          with.lambda = TRUE,
          init.attempt = init.attempt
        )
      }
      ## TODO: hydromad.getOption("catch.errors")
      pcheck <- try(tfParsCheck(pars))
      if (!isTRUE(pcheck)) {
        return(pcheck)
      }
    }
    oldU <- NULL
    oldpars <- pars
    i <- 1
    repeat {
      if (init.U && (i == 1)) {
        ## we already have an estimate for U
      } else {
        U <- lambda.inverse.sim(Q,
          P = P,
          pars = pars, delay = delay,
          use.Qm = use.Qm,
          rises.only = rises.only
        )
      }
      ## inverse.sim should be in sync with P not Q,
      ## so delay still applies
      fnName <- paste("lambda", fit.method, "fit", sep = ".")
      modCall <- quote(FUN(cbind(Q = Q, U = U),
        delay = delay, ...
      ))
      # lambda.init = pars[["lambda"]], ...))
      modCall[[1]] <- as.symbol(fnName)
      mod <- eval(modCall)
      if (!inherits(mod, "tf")) {
        return(mod)
      }

      pars <- coef(mod)

      if (trace) {
        message(paste(" iteration = ", i))
        print(pars)
      }

      if (!is.null(oldU)) {
        delta <- sum(abs(U - oldU), na.rm = TRUE) / sumQ
        if (delta < rel.tolerance) break
      }

      if (i >= max.iterations) {
        warning("reached maximum number of iterations")
        break
      }
      oldU <- U
      oldpars <- pars
      i <- i + 1
    }

    mod$call <- match.call()
    mod
  }
