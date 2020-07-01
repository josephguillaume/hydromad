## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

expuh.ls.fit <-
  function(DATA,
           order = hydromad.getOption("order"),
           delay = hydromad.getOption("delay"),
           quiet = FALSE,
           ...) {
    model <- armax.ls.fit(DATA, order = order, delay = delay, ...)
    if (!inherits(model, "tf")) {
      return(model)
    }
    n <- order[1]
    poles <- arToPoles(coef(model)[seq_len(n)])
    eps <- sqrt(.Machine$double.eps)
    badpoles <- (Re(poles) < -eps) | (abs(Im(poles)) > eps)
    if (any(badpoles)) {
      if (!quiet) {
        message("armax fitted poles are non-physical; re-fitting with constraints")
      }
      model <-
        fitWithPoleConstraints(DATA,
          fitfun = armax.ls.fit, poles = poles,
          order = order, delay = delay, ...
        )
    }
    if (!inherits(model, "tf")) {
      return(model)
    }
    model$coefficients <- coef(model, "tau,v")
    model
  }


expuh.sriv.fit <-
  function(DATA,
           order = hydromad.getOption("order"),
           delay = hydromad.getOption("delay"),
           quiet = FALSE,
           ...) {
    model <- armax.sriv.fit(DATA, order = order, delay = delay, ...)
    if (!inherits(model, "tf")) {
      return(model)
    }
    n <- order[1]
    poles <- arToPoles(coef(model)[seq_len(n)])
    eps <- sqrt(.Machine$double.eps)
    badpoles <- (Re(poles) < -eps) | (abs(Im(poles)) > eps)
    if (any(badpoles)) {
      if (!quiet) {
        message("armax fitted poles are non-physical; re-fitting with constraints")
      }
      model <-
        fitWithPoleConstraints(DATA,
          fitfun = armax.sriv.fit, poles = poles,
          order = order, delay = delay, ...
        )
    }
    if (!inherits(model, "tf")) {
      return(model)
    }
    model$coefficients <- coef(model, "tau,v")
    model
  }


expuh.inverse.fit <-
  function(DATA,
           order = hydromad.getOption("order"),
           delay = hydromad.getOption("delay"),
           ...) {
    ## TODO: can do this directly?
    model <- armax.inverse.fit(DATA, order = order, delay = delay, ...)
    if (!inherits(model, "tf")) {
      return(model)
    }
    # model$vcov <- vcov(model)
    model$coefficients <- coef(model, "tau,v")
    model
  }

fitWithPoleConstraints <-
  function(DATA, fitfun, poles, ...,
           control = as.list(hydromad.getOption("optim.control.expuh"))) {
    ## non-physical roots; constrain them and fit again
    poles <- abs(Re(poles))
    poles <- pmax(poles, 1e-5)
    logpol0 <- log(poles)
    model <- structure("failed to fit expuh routing with constrained roots",
      class = "try-error"
    )
    bestVal <- Inf
    optFun <- function(logpol, ...) {
      pol <- exp(logpol)
      if (isTRUE(hydromad.getOption("catch.errors"))) {
        thisMod <- try(fitfun(DATA, ..., fixed.ar = polesToAr(pol)))
      } else {
        thisMod <- fitfun(DATA, ..., fixed.ar = polesToAr(pol))
      }
      if (!isValidModel(thisMod)) {
        return(NA)
      }
      val <- sum(abs(residuals(thisMod)), na.rm = TRUE)
      if (val < bestVal) {
        model <<- thisMod
      }
      val
    }
    if (!isTRUE(hydromad.getOption("catch.errors.optim"))) {
      try <- force
    } ## i.e. skip the try()
    # ans <- try(nlminb(logpol0, optFun, control = control, ...))
    ans <- try(optim(logpol0, optFun, control = control, ...))
    if (inherits(ans, "try-error")) {
      return(ans)
    }
    if (ans$convergence != 0) {
      msg <- paste("While re-fitting non-physical poles:", toString(ans$msg))
      if (!isTRUE(hydromad.getOption("quiet"))) {
        warning(msg)
      }
      model$msg <- msg
    }
    model
  }
