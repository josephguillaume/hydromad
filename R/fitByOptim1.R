## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


fitByOptim1 <-
  function(MODEL,
           objective = hydromad.getOption("objective"),
           tol = .Machine$double.eps^0.25) {
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
    if (length(parlist) > 1) {
      stop("this function can only handle a single free parameter")
    }
    lower <- sapply(parlist, min)
    upper <- sapply(parlist, max)
    bestModel <- MODEL
    bestFunVal <- -Inf
    objseq <- rep(NA_real_, 100)
    i <- 0
    do_optim1 <- function(pars) {
      i <<- i + 1
      names(pars) <- names(parlist)
      thisMod <- update(MODEL, newpars = pars)
      if (!isValidModel(thisMod)) {
        return(NA)
      }
      thisVal <- objFunVal(thisMod, objective = objective)
      objseq[i] <<- thisVal
      if (isTRUE(thisVal > bestFunVal)) {
        bestModel <<- thisMod
        bestFunVal <<- thisVal
      }
      ## optimize does minimisation, so:
      return(-thisVal)
    }
    ans <- optimize(do_optim1,
      lower = lower, upper = upper,
      tol = tol
    )
    bestModel$funevals <- i
    bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
    bestModel$objective <- objective
    bestModel$objseq <- objseq[1:i]
    bestModel$fit.call <- match.call()
    bestModel$fit.result <- ans
    return(bestModel)
  }
