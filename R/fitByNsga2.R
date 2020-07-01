fitByNsga2 <- function(MODEL, objective = hydromad.getOption("objective"),
                       control = hydromad.getOption("nsga2.control")) {
  if (!requireNamespace("mco")) stop("package mco is required for fitByNsga2")
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
  #     if (!isTRUE(hydromad.getOption("trace")))
  #         control$trace <- FALSE
  lower <- sapply(parlist, min)
  upper <- sapply(parlist, max)
  bestModel <- MODEL
  bestFunVal <- -Inf
  do_nsga2 <- function(pars) {
    names(pars) <- names(parlist)
    thisMod <- update(MODEL, newpars = pars)
    if (!isValidModel(thisMod)) {
      return(1e+08)
    }
    thisVal <- objFunVal(thisMod, objective = objective)
    if (isTRUE(thisVal > bestFunVal)) {
      bestModel <<- thisMod
      bestFunVal <<- thisVal
    }
    return(-thisVal)
  }
  args <- modifyList(control, list(
    fn = do_nsga2, idim = length(parlist), odim = 1,
    lower.bounds = lower, upper.bounds = upper
  ))
  ans <- do.call(mco::nsga2, args)
  bestModel$funevals <- NA ## TODO
  bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
  bestModel$objective <- objective
  bestModel$fit.call <- match.call()
  bestModel$fit.result <- ans
  return(bestModel)
}
