## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


## NOTE - currently all the tf() stuff converts data to ts

tf <-
  function(DATA = zoo(),
           pars = numeric(),
           delay = 0, # hydromad.getOption("delay"),
           warmup = hydromad.getOption("warmup"),
           initX = TRUE,
           na.action = na.pass,
           epsilon = TRUE) {
    if (is.na(delay)) delay <- 0
    if (missing(DATA)) {
      ## allow to define tf without data, for later predict()ion
      DATA <- as.ts(cbind(
        U = rep(0, length = warmup + 1),
        Q = rep(0, length = warmup + 1)
      ))
    }
    ## get data into the right form
    DATA <- as.ts(DATA)
    DATA <- na.action(DATA)
    stopifnot(c("U", "Q") %in% colnames(DATA))
    ## the model object
    obj <- list(call = match.call())
    obj$data <- DATA
    class(obj) <- c("tf")
    obj <- update(obj,
      pars = pars, delay = delay, warmup = warmup,
      initX = initX, epsilon = epsilon
    )
    if (!inherits(obj, "tf")) {
      return(obj)
    }
    obj$call <- match.call()
    return(obj)
  }

update.tf <-
  function(object, newdata = NULL, pars, ...) {
    upcall <- match.call()
    nm <- names(upcall)
    if (!is.null(nm)) {
      nm <- nm[!(nm %in% c("", "object", "newdata"))]
      if (length(nm) > 0) {
        object$call[nm] <- upcall[nm]
      }
    }
    ## update arguments (delay, warmup, initX, epsilon)
    dots <- list(...)
    stopifnot(all(names(dots) %in%
      c("delay", "warmup", "initX", "epsilon")))
    object <- modifyList(object, dots)
    ## update parameters
    if (!missing(pars)) {
      if (length(pars) == 0) {
        object$order <- c(n = 0, m = 0)
        object$coefficients <- numeric()
      } else {
        pars <- tfParsConvert(pars, "a,b")
        ## check parameters
        if (isTRUE(hydromad.getOption("catch.errors"))) {
          pcheck <- try(tfParsCheck(pars),
            silent = !hydromad.getOption("trace")
          )
          if (!isTRUE(pcheck)) {
            return(pcheck)
          }
        } else {
          tfParsCheck(pars)
        }
        ## model order
        parSymbols <- gsub("_.*", "", names(pars))
        a <- pars[parSymbols == "a"]
        b <- pars[parSymbols == "b"]
        a <- stripzeros(a)
        b <- stripzeros(b, up.to = 1)
        n <- length(a)
        m <- length(b) - 1
        object$order <- c(n = n, m = m)
        extraPars <- pars[(parSymbols %in% c("a", "b")) == FALSE]
        object$coefficients <- c(a, b, extraPars)
      }
    }
    ## update DATA
    if (!is.null(newdata)) {
      newdata <- as.ts(newdata)
      stopifnot("Q" %in% colnames(newdata))
      object$data <- newdata
    }

    initX <- if (!is.null(object$initX)) object$initX else 0
    if (length(coef(object)) > 0) {
      ## take starting value of filter from data
      if (isTRUE(initX)) {
        Q <- object$data[, "Q"]
        if (any(is.finite(Q[1:10]))) {
          initX <- min(Q[1:10], na.rm = TRUE)
        }
      }
      ## determine data resolution from data
      epsilon <- object$epsilon
      if (isTRUE(epsilon)) {
        epsilon <- hydromad.getOption("sim.epsilon")
        Q <- object$data[, "Q"]
        if (any(Q < epsilon, na.rm = TRUE)) {
          isNonZero <- is.finite(Q) & (Q > 0)
          if (any(isNonZero)) {
            eps.data <- min(Q[isNonZero])
            if (eps.data < 1) {
              epsilon <- eps.data / 10
            }
          }
        }
        object$epsilon <- epsilon
      }
    }
    ## run model
    X <- predict(object, init = initX, epsilon = object$epsilon)
    object$fitted.values <- X
    return(object)
  }

predict.tf <-
  function(object,
           newdata = NULL,
           ...) {
    if (is.null(newdata)) {
      newdata <- object$data
    }
    ## get data into the right form
    newdata <- as.ts(newdata)
    if (NCOL(newdata) > 1) {
      newdata <- newdata[, "U"]
    }
    ## simulate
    armax.sim(newdata,
      pars = coef(object),
      delay = object$delay,
      ...
    )
  }

print.tf <-
  function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nUnit Hydrograph / Linear Transfer Function\n")
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    n <- x$order["n"]
    m <- x$order["m"]
    cat("Order: (n=", n, ", m=", m, ")  Delay: ", x$delay, "\n", sep = "")
    cat("ARMAX Parameters:\n")
    print.default(coef(x), digits = digits, print.gap = 2)
    cat("Exponential component parameters:\n")
    pars.tv <- try(coef(x, "tau,v"), silent = TRUE)
    print.default(pars.tv, digits = digits, print.gap = 2)
    cat("TF Structure: ")
    cat(describeTF(coef(x)))
    cat("\n\n")
    invisible(x)
  }

coef.tf <- function(object, form = c("a,b", "tau,v", "alpha,beta"), ...) {
  coef <- object$coefficients
  if (!missing(form)) {
    form <- match.arg(form)
    return(tfParsConvert(coef, form))
  }
  coef
}

xyplot.tf <- function(x, data = NULL, ...) {
  xyplot.hydromad(x, ...)
}

ssg.tf <-
  function(object, ...) {
    theta <- coef(object)
    #    if (length(theta) == 0) {
    #        ## estimate scale by mass balance
    #        Q <- object$data[,"Q"]
    #        U <- object$data[,"U"]
    #        ok <- complete.cases(Q, U)
    #        return(sum(U[ok]) / sum(Q[ok]))
    #    } else {
    ssg.tf.coef(theta)
    #    }
  }


normalise.tf <-
  function(object, ...) {
    ## TODO: do not run model again, just scale output?
    update(object, pars = normalise.tf.coef(coef(object)))
  }

fitted.tf <- function(object, ..., all = FALSE) {
  tmp <- object$fitted.values
  if (length(tmp) == 0) {
    return(tmp)
  }
  if (all) {
    return(tmp)
  }
  stripWarmup(tmp, object$warmup)
}

residuals.tf <- function(object, ..., all = FALSE) {
  f <- fitted(object, all = TRUE)
  if (length(f) == 0) {
    return(f)
  }
  tmp <- (object$data[, "Q"] - f)
  if (all) {
    return(tmp)
  }
  stripWarmup(tmp, object$warmup)
}

observed.tf <- function(object, ..., all = FALSE) {
  ## observed.default will work, but this may be slightly faster
  tmp <- object$data[, "Q"]
  if (all) {
    return(tmp)
  }
  stripWarmup(tmp, object$warmup)
}

vcov.tf <- function(object, ...) {
  object$cov.mat
} ## may be NULL

deviance.tf <- stats:::deviance.lm

isValidModel.tf <- function(object, ...) {
  (inherits(object, "tf") &&
    (isTRUE(try(tfParsCheck(coef(object)), silent = TRUE))))
}
