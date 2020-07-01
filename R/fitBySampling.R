## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


fitBySampling <-
  function(MODEL,
           objective = hydromad.getOption("objective"),
           samples = hydromad.getOption("fit.samples"),
           sampletype = c("latin.hypercube", "random", "all.combinations")) {
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
    ## generate parameter sets
    psets <- parameterSets(parlist, samples = samples, method = sampletype)
    bestModel <- MODEL
    bestFunVal <- -Inf
    objseq <- rep(NA_real_, NROW(psets))
    for (i in seq(NROW(psets))) {
      thisPars <- as.list(psets[i, , drop = FALSE])
      if (isTRUE(hydromad.getOption("trace"))) {
        run_name <- paste(names(thisPars), signif(unlist(thisPars), 3),
          sep = "=", collapse = ", "
        )
        message(run_name)
      }
      thisMod <- update(MODEL, newpars = thisPars)
      if (!isValidModel(thisMod)) {
        next
      }
      thisVal <- objFunVal(thisMod,
        objective = objective,
        nan.ok = hydromad.getOption("catch.errors")
      )
      objseq[i] <- thisVal
      if (isTRUE(thisVal > bestFunVal)) {
        bestModel <- thisMod
        bestFunVal <- thisVal
      }
    }
    bestModel$funevals <- NROW(psets)
    bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
    bestModel$objective <- objective
    bestModel$fit.call <- match.call()
    bestModel$fit.result <- list(objseq = objseq, psets = psets)
    bestModel
  }
