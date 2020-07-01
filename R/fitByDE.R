## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


fitByDE <-
  function(MODEL,
           objective = hydromad.getOption("objective"),
           control = hydromad.getOption("de.control")) {
    if (!requireNamespace("DEoptim")) stop("package DEoptim is required for fitByDE")
    control <- do.call(DEoptim::DEoptim.control, control)
    start_time <- proc.time()
    objective <- buildCachedObjectiveFun(objective, MODEL)
    parlist <- as.list(coef(MODEL, warn = FALSE))
    ## remove any missing parameters
    isok <- sapply(parlist, function(x) !any(is.na(x)))
    parlist <- parlist[isok]
    ## check which parameters are uniquely specified
    isfixed <- (sapply(parlist, length) == 1)
    if (all(isfixed)) {
      warning("all parameters are fixed, so can not fit")
      return(MODEL)
    }
    ## remove any fixed parameters
    parlist <- parlist[!isfixed]
    if (!isTRUE(hydromad.getOption("trace"))) {
      control$trace <- FALSE
    }
    lower <- sapply(parlist, min)
    upper <- sapply(parlist, max)
    bestModel <- MODEL
    bestFunVal <- -Inf
    do_de <- function(pars) {
      names(pars) <- names(parlist)
      thisMod <- update(MODEL, newpars = pars)
      if (!isValidModel(thisMod)) {
        return(1e8)
      }
      thisVal <- objFunVal(thisMod, objective = objective)
      if (isTRUE(thisVal > bestFunVal)) {
        bestModel <<- thisMod
        bestFunVal <<- thisVal
      }
      ## DEoptim does minimisation, so:
      return(-thisVal)
    }
    ans <- DEoptim::DEoptim(do_de,
      lower = lower, upper = upper,
      control = control
    )
    bestModel$funevals <- ans$optim$nfeval
    bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
    bestModel$objective <- objective
    bestModel$fit.call <- match.call()
    bestModel$fit.result <- ans
    return(bestModel)
  }
