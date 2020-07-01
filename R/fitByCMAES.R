fitByCMAES <- function(MODEL, objective = hydromad.getOption("objective"),
                       control = hydromad.getOption("cmaes.control"), vcov = FALSE) {
  if (!requireNamespace("cmaes")) stop("package cmaes is required for fitByCMAES")
  start_time <- proc.time()
  objective <- buildCachedObjectiveFun(objective, MODEL)
  parlist <- as.list(coef(MODEL, warn = FALSE))
  isok <- sapply(parlist, function(x) !any(is.na(x)))
  parlist <- parlist[isok]
  isfixed <- (sapply(parlist, length) == 1)
  if (all(isfixed)) {
    warning("all parameters are fixed, so can not fit")
    return(MODEL)
  }
  parlist <- parlist[!isfixed]
  ## TODO: if (isTRUE(hydromad.getOption("trace")))
  lower <- sapply(parlist, min)
  upper <- sapply(parlist, max)
  initpars <- sapply(parlist, mean)
  bestModel <- MODEL
  bestFunVal <- -Inf
  do_cmaes <- function(pars) {
    thisMod <- update(MODEL, newpars = pars)
    if (!isValidModel(thisMod)) {
      return(NA)
    }
    thisVal <- objFunVal(thisMod, objective = objective)
    if (isTRUE(thisVal > bestFunVal)) {
      bestModel <<- thisMod
      bestFunVal <<- thisVal
    }
    return(-thisVal)
  }
  ans <- cmaes::cma_es(initpars, do_cmaes,
    lower = lower, upper = upper,
    control = control
  )
  if (ans$convergence != 0) {
    if (!isTRUE(hydromad.getOption("quiet"))) {
      warning(ans$message)
    }
    bestModel$msg <- ans$message
  }
  bestModel$funevals <- ans$counts[1]
  bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
  bestModel$objective <- objective
  if (vcov) {
    warning("vcov not yet implemented")
  }
  bestModel$fit.call <- match.call()
  bestModel$fit.result <- ans
  return(bestModel)
}
