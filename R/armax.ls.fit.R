## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##




#' Estimate transfer function models by Least Squares.
#'
#' Calibrate unit hydrograph transfer function models (\code{\link{armax}} or
#' \code{\link{expuh}}) using Least Squares with prefiltering.
#'
#' In normal usage, one would not call these functions directly, but rather
#' specify the routing fitting method for a \code{\link{hydromad}} model using
#' that function's \code{rfit} argument. E.g. to specify fitting an
#' \code{expuh} routing model by least squares one could write
#'
#' \code{hydromad(..., routing = "expuh", rfit = "ls")}
#'
#' which uses the default order, \code{hydromad.getOption("order")}, or
#'
#' \code{hydromad(..., routing = "expuh", rfit = list("ls", order = c(2,1)))}.
#'
#' @importFrom stats na.pass na.omit ts.intersect lm.wfit lm.fit residuals var
#'
#' @name armax.ls.fit
#' @param DATA a \code{\link{ts}}-like object with named columns: \describe{
#' \item{list("U")}{ observed input time series. } \item{list("Q")}{ observed
#' output time series. } }
#' @param order the transfer function order. See \code{\link{armax}}.
#' @param delay delay (lag time / dead time) in number of time steps. If
#' missing, this will be estimated from the cross correlation function.
#' @param prefilter placeholder
#' @param warmup placeholder
#' @param normalise placeholder
#' @param fixed.ar placeholder
#' @param weights placeholder
#' @param initX placeholder
#' @param na.action placeholder
#' @param trace placeholder
#' (i.e. negative or imaginary poles) are detected.
#' @return a \code{tf} object, which is a list with components
#' \item{coefficients}{ the fitted parameter values.} \item{fitted.values}{ the
#' fitted values. } \item{residuals}{ the residuals. } \item{delay}{ the
#' (possibly fitted) delay time. }
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{armax}}, \code{\link{expuh}},
#' \code{\link{armax.sriv.fit}}, \code{\link{arima}}
#' @references Jakeman
#' @keywords ts
#' @examples
#'
#' U <- ts(c(0, 0, 0, 1, rep(0, 30), 1, rep(0, 20)))
#' Y <- expuh.sim(lag(U, -1), tau_s = 10, tau_q = 2, v_s = 0.5, v_3 = 0.1)
#' set.seed(0)
#' Yh <- Y * rnorm(Y, mean = 1, sd = 0.2)
#' fit1 <- armax.ls.fit(ts.union(U = U, Q = Yh),
#'   order = c(2, 2), warmup = 0
#' )
#' fit1
#' xyplot(ts.union(observed = Yh, fitted = fitted(fit1)),
#'   superpose = TRUE
#' )
#' @export armax.ls.fit
armax.ls.fit <-
  function(DATA,
           order = hydromad.getOption("order"),
           delay = hydromad.getOption("delay"),
           prefilter = hydromad.getOption("prefilter"),
           warmup = hydromad.getOption("warmup"),
           normalise = FALSE,
           fixed.ar = NULL,
           weights = NULL,
           initX = TRUE,
           na.action = na.pass,
           trace = hydromad.getOption("trace")) {
    ## get data into the right form
    DATA <- as.ts(DATA)
    DATA <- na.action(DATA)
    stopifnot(c("U", "Q") %in% colnames(DATA))
    if (is.na(delay)) {
      delay <- estimateDelay(DATA, plot = FALSE)
    }
    ## check values
    stopifnot(length(order) == 2)
    if (is.null(names(order))) names(order) <- c("n", "m")
    n <- order[["n"]]
    m <- order[["m"]]
    warmup0 <- warmup
    warmup <- max(n, m + delay, warmup) ## used in fitting only
    stopifnot(all(is.finite(fixed.ar)))

    DATAf <- DATA
    if (!is.null(weights)) {
      weights <- as.ts(weights)
      weights <- na.action(weights)
      ## include weights in data matrix, it will be filtered
      nm <- colnames(DATA)
      DATAf <- ts.intersect(DATAf, weights)
      colnames(DATAf) <- c(nm, "weights")
    }
    # if (!identical(prefilter, FALSE)) {
    #    nm <- colnames(DATAf)
    #    DATAf <- filter(DATAf, filter = prefilter, method = "recursive")
    #    colnames(DATAf) <- nm
    # }
    ## construct model matrix
    if (n > 0) {
      z <- ts.intersect(
        Y = DATAf[, "Q"],
        Q = lagbind(DATAf[, "Q"], lags = -(1:n), all = FALSE),
        U = lagbind(DATAf[, "U"], lags = -(0:m) - delay, all = FALSE)
      )
    } else {
      z <- ts.intersect(
        Y = DATAf[, "Q"],
        U = lagbind(DATAf[, "U"], lags = -(0:m) - delay, all = FALSE)
      )
    }
    ## align weights
    if (!is.null(weights)) {
      w <- window(DATAf[, "weights"], start = start(z), end = end(z))
    }

    ## strip off warmup period; this also converts to matrix
    # z <- z[-(1:warmup), ]
    z <- as.matrix(z)

    ## remove any remaining NAs from data matrix
    z[!is.finite(z)] <- NA
    z <- z[complete.cases(z), ]

    ## HEURISTIC METHOD FOR PREFILTER:
    if (identical(prefilter, FALSE)) {
      prefilter <- list(0)
    } else {
      if (!is.list(prefilter)) {
        prefilter <- list(prefilter)
      }
    }

    obj <- NULL
    fit.error <- Inf
    which.pf <- NA

    for (i in seq_along(prefilter)) {
      ## apply i'th candidate prefilter
      pf <- prefilter[[i]]
      ## handle prefilter-generating function
      if (is.function(pf)) {
        pf <- pf(DATA)
      }
      z2 <- filter_ok(z, filter = pf, method = "recursive")
      ## set column names for parameters
      a_names <- if (n > 0) paste("a", 1:n, sep = "_")
      b_names <- paste("b", 0:m, sep = "_")
      colnames(z2) <- c("y", a_names, b_names)
      ## apply fixed parameters if any
      if (length(fixed.ar) > 0) {
        ## substract fixed predictors
        z2[, 1] <- z2[, 1] - z2[, 1 + (1:n)] %*% fixed.ar
        ## and remove them from model matrix
        z2 <- z2[, -(1 + (1:n)), drop = FALSE]
      }
      z2 <- z2[-(1:warmup), , drop = FALSE]
      ## fit the model
      if (!is.null(weights)) {
        w2 <- filter_ok(w, filter = pf, method = "recursive")
        w2 <- w2[-(1:warmup)]
        fit.f <- lm.wfit(x = z2[, -1, drop = FALSE], y = z2[, 1], w = w2)
      } else {
        fit.f <- lm.fit(x = z2[, -1, drop = FALSE], y = z2[, 1])
      }
      ## extract calibrated parameters and re-name them
      theta <- coef(fit.f)
      theta[!is.finite(theta)] <- 0 ## obscure case: all(Q*U==0)
      if (length(fixed.ar) > 0) {
        names(fixed.ar) <- a_names
        theta <- c(fixed.ar, theta)
      }
      if (n > 0) {
        theta[1:n] <- stabiliseAR(theta[1:n])
      }
      if (normalise) {
        theta <- normalise.tf.coef(theta)
      }
      theta[!is.finite(theta)] <- 0
      ## create fitted model object (includes simulation)
      ## this will check that solution is OK
      obji <- tf(DATA, pars = theta, delay = delay, warmup = warmup0, initX = initX)
      if (!inherits(obji, "tf")) {
        next
      }
      errori <- mean(abs(residuals(obji)), na.rm = TRUE)
      if (errori < fit.error) {
        obj <- obji
        fit.error <- errori
        which.pf <- i
      }
    }
    if (!inherits(obj, "tf")) {
      return(obji)
    }
    if (trace && (length(prefilter) > 1)) {
      message("chosen prefilter: ", which.pf)
    }
    p <- fit.f$rank
    cov.unscaled <- try(chol2inv(fit.f$qr$qr[1:p, 1:p, drop = FALSE]))
    obj$call <- match.call()
    obj$weights <- weights
    obj$prefilter <- prefilter[[which.pf]]
    ## parameter covariance matrix
    obj$sigma2 <- var(residuals(obj), na.rm = TRUE)
    if (!inherits(cov.unscaled, "try-error")) {
      cov.mat <- obj$sigma2 * cov.unscaled
      colnames(cov.mat) <- rownames(cov.mat) <-
        names(na.omit(coef(fit.f)))
      obj$cov.mat <- cov.mat
    }
    obj
  }
