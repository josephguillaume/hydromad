## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


armax.sriv.fit <-
  function(DATA,
           order = hydromad.getOption("order"),
           delay = hydromad.getOption("delay"),
           noise.order = hydromad.getOption("riv.noise.order"),
           fixed.ar = NULL,
           ...,
           fallback = TRUE,
           na.action = na.pass,
           epsilon = hydromad.getOption("sriv.epsilon"),
           max.iterations = hydromad.getOption("sriv.iterations")) {
    ## get data into the right form
    DATA <- as.ts(DATA)
    DATA <- na.action(DATA)
    stopifnot(c("U", "Q") %in% colnames(DATA))
    if (is.na(delay)) {
      delay <- estimateDelay(DATA, plot = FALSE)
    }

    ## first fit by least squares
    init.model <-
      armax.ls.fit(DATA,
        order = order, delay = delay,
        fixed.ar = fixed.ar,
        ...
      )
    if (!inherits(init.model, "tf")) {
      return(init.model)
    }

    obj <- do_srivfit(DATA,
      init.model = init.model,
      order = order,
      noise.order = noise.order,
      fixed.ar = fixed.ar,
      ...,
      fallback = fallback,
      epsilon = epsilon,
      max.iterations = max.iterations
    )
    if (inherits(obj, "tf")) {
      obj$call <- match.call()
    }
    obj
  }

do_srivfit <-
  function(DATA,
           init.model,
           prefilter = hydromad.getOption("prefilter"),
           order = hydromad.getOption("order"),
           noise.order = hydromad.getOption("riv.noise.order"),
           fixed.ar = NULL,
           warmup = hydromad.getOption("warmup"),
           normalise = FALSE,
           initX = TRUE,
           na.action = na.pass,
           fallback = TRUE,
           epsilon = hydromad.getOption("sriv.epsilon"),
           max.iterations = hydromad.getOption("sriv.iterations"),
           trace = hydromad.getOption("trace")) {
    stopifnot(c("U", "Q") %in% colnames(DATA))
    ## check values
    stopifnot(inherits(init.model, "tf"))
    delay <- init.model$delay
    stopifnot(length(order) == 2)
    if (is.null(names(order))) names(order) <- c("n", "m")
    n <- order[["n"]]
    m <- order[["m"]]
    # if (n == 0) {
    #    warning("SRIV algorithm only works with models with an AR component (n>0)")
    #    return(init.model)
    # }
    stopifnot(epsilon > 0)
    warmup0 <- warmup
    warmup <- max(n, m + delay, warmup) ## used in fitting only
    ## observed variance (excludes warmup period)
    if (trace) obs.var <- var(DATA[-(1:warmup), "Q"], na.rm = TRUE)
    theta <- coef(init.model)
    initorder <- init.model$order
    if (any(init.model$order < order)) {
      ## this can happen if e.g. one AR term is (near)zero:
      ## the returned model from armax.ls.fit has reduced order.
      ## pad out parameter vector with zeros.
      initn <- init.model$order[["n"]]
      initm <- init.model$order[["m"]]
      theta <- rep(0, n + m + 1)
      theta[seq_len(initn)] <- coef(init.model)[seq_len(initn)]
      theta[n + seq_len(initm)] <- coef(init.model)[initn + seq_len(initm)]
    }
    theta.prev <- theta
    a.hat <- theta[1:n]
    if (n == 0) {
      a.hat <- 0
    }
    mask <- rep(TRUE, length = length(theta))
    if (length(fixed.ar) > 0) {
      mask[1:n] <- FALSE
    }
    X <- fitted(init.model, all = TRUE)
    obj <- init.model
    ## working data matrix
    DATA <- DATA[, c("U", "Q", "U")] ## last column will store X
    converged <- FALSE
    for (iteration in seq(max.iterations)) {
      ## SRIV parameter estimate...
      DATAf <- DATA
      stopifnot(tsp(X) == tsp(DATA))
      DATAf[, 3] <- X ## NOTE X[] vector has already been aligned with data
      ## DATAf <- ts.intersect(DATA[,c("U","Q")], X=X) ## SAFER
      if (!identical(prefilter, FALSE)) {
        if (any(noise.order > 0)) {
          ## RIV
          nn <- noise.order[1]
          nm <- noise.order[2]
          i.arma <- arima(X - DATA[, "Q"],
            order = c(nn, 0, nm),
            include.mean = FALSE
          )
          nar <- coef(i.arma)[seq_len(nn)]
          nma <- coef(i.arma)[nn + seq_len(nm)]
          ## apply model AR and inverse of noise
          ## TODO: is this correct?
          DATAf <- filter_ok(DATAf, filter = nar, sides = 1)
          DATAf <- filter_ok(DATAf, filter = nma, method = "recursive")
          DATAf <- filter_ok(DATAf, filter = a.hat, method = "recursive")
        } else {
          ## SRIV: auto-regressive filter only
          DATAf <- filter_ok(DATAf, filter = a.hat, method = "recursive")
        }
      }
      colnames(DATAf) <- c("U", "Q", "X")

      if (hydromad.getOption("pure.R.code") == FALSE) {
        ## TODO... handle NA / NaN / Inf
        ans <- .C(sriv_system,
          as.double(DATAf[, "U"]),
          as.double(DATAf[, "Q"]),
          as.double(DATAf[, "X"]),
          as.integer(NROW(DATAf)),
          as.integer(warmup),
          as.integer(order),
          as.integer(delay),
          xz = double((n + m + 1)^2),
          xy = double((n + m + 1)),
          xx = double((n + m + 1)^2),
          NAOK = TRUE, DUP = FALSE, PACKAGE = "hydromad"
        )
        xz <- matrix(ans$xz, ncol = (n + m + 1))
        xx <- matrix(ans$xx, ncol = (n + m + 1))
        xQ <- ans$xy
      } else {
        ## implementation in R for cross-checking
        Uf <- DATAf[, "U"]
        Qf <- DATAf[, "Q"]
        Xf <- DATAf[, "X"]
        ## z: regressors { Qf[k-1] ... Qf[k-n], Uf[k] ... Uf[k-m] }
        ## x: instrumental variable { Xf[k-1] ... Xf[k-n], Uf[k] ... Uf[k-m] }
        ## form the system as a matrix outer product (x z')
        ## row t of z is the regressor vector z at time t
        z <- ts.intersect(
          Q = lagbind(Qf, lags = -(1:n), all = FALSE),
          U = lagbind(Uf, lags = -(0:m) - delay, all = FALSE)
        )
        ## row t of x is the instrumental variable x at time t
        x <- ts.intersect(
          X = lagbind(Xf, lags = -(1:n), all = FALSE),
          U = z[, -(1:n)]
        )
        ## strip off warmup period
        x <- stripWarmup(x, warmup)
        z <- stripWarmup(z, warmup)
        # x <- window(x, start=t_warm, end=t_end)
        # z <- window(z, start=t_warm, end=t_end)
        ## each row of x_and_z is c(x_t, z_t) for time t
        x_and_z <- ts.intersect(x = x, z = z)
        x_and_z <- x_and_z[complete.cases(x_and_z), ]
        x.i <- 1:ncol(x)
        z.i <- 1:ncol(z) + ncol(x)
        ## compute outer products of x_t and z_t for all times t
        xz <- apply(x_and_z, ROWS <- 1, function(data_t) {
          data_t[x.i] %*% t(data_t[z.i])
        })
        ## (now each column of xz holds the entries of the outer product)
        ## sum over time (rowSums) and form back into a matrix
        xz <- matrix(rowSums(xz), nrow = (n + m + 1))
        colnames(xz) <- colnames(z)
        rownames(xz) <- colnames(x)
        ## form the system "response" as product (x Q)
        ## note ts times are aligned automatically in `*`
        xQ <- colSums(x * Qf, na.rm = TRUE)
        ## form the information matrix as (x x')
        xx <- apply(x, ROWS <- 1, function(x_t) {
          x_t %*% t(x_t)
        })
        xx <- matrix(rowSums(xx), nrow = (n + m + 1))
        rownames(xx) <- colnames(xx) <- colnames(x)
      }
      ## solve the system for parameter set { a_1 ... a_n, b_0 ... b_m }
      ##
      ## apply fixed parameters if any
      if (length(fixed.ar) > 0) {
        ## substract fixed predictors
        xQ <- drop(xQ - xz[, (1:n)] %*% fixed.ar)
        ## and remove them from model matrix
        xz <- xz[, -(1:n), drop = FALSE]
        theta <- try(qr.solve(xz, xQ),
          silent = !hydromad.getOption("trace")
        )
      } else {
        theta <- try(solve(xz, xQ),
          silent = !hydromad.getOption("trace")
        )
      }
      if (inherits(theta, "try-error")) {
        if (fallback) {
          warning(theta)
          obj$warning <- theta
          theta <- theta.prev
          break
        } else {
          return(theta)
        }
      }
      if (length(fixed.ar) > 0) {
        theta <- c(fixed.ar, theta)
      }
      names(theta) <- names(theta.prev)
      if (n > 0) {
        theta[1:n] <- stabiliseAR(theta[1:n])
        a.hat <- theta[1:n]
      }
      if (trace) {
        if (iteration == 1) {
          print(round(theta.prev, 8))
        }
        print(round(theta, 8))
      }
      ## create fitted model object (includes simulation)
      ## this will check that solution is OK
      obj <- tf(DATA,
        pars = theta, delay = delay,
        warmup = warmup0, initX = initX
      )
      if (!inherits(obj, "tf")) {
        return(obj)
      }
      X <- fitted(obj, all = TRUE)
      ## stop if all parameters converged within epsilon
      deltas <- abs((theta - theta.prev) / theta)
      if (mean(deltas[mask], na.rm = TRUE) < epsilon) {
        converged <- TRUE
        break
      }
      ## go on
      theta.prev <- theta
    }
    if (!converged) {
      warning("did not converge after ", iteration, " iterations")
    }
    resid.var <- var(residuals(obj), na.rm = TRUE)
    info.mat <- xx / resid.var
    cov.mat <- matrix()
    cm.solve <- try(
      {
        cov.mat <- solve(info.mat) ## symmetric -- use chol?
      },
      silent = !hydromad.getOption("trace")
    )
    if (!inherits(cm.solve, "try-error")) colnames(cov.mat) <- rownames(cov.mat) <- names(theta)
    if (trace) {
      if (converged) {
        message("converged after ", iteration, " iterations")
      }
      ## report performance stats
      arpe <- mean(diag(cov.mat) / (theta^2))
      message("ARPE (%): ", format(arpe * 100))
      message(
        "1 - (resid.var/obs.var):  ",
        format(1 - (resid.var / obs.var))
      )
    }
    obj$converged <- converged
    obj$iteration <- iteration
    obj$sigma2 <- resid.var
    obj$cov.mat <- cov.mat
    if (normalise) obj <- normalise.tf(obj)
    obj
  }
