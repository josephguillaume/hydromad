## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

fitByDream <-
  function(MODEL,
           loglik = hydromad.getOption("loglik"),
           control = hydromad.getOption("dream.control"),
           vcov = TRUE, save = NULL) {
    if (!requireNamespace("dream")) stop('package dream is required for fitByDream.\n  Use: install.packages("dream", repos="http://hydromad.catchment.org")')
    start_time <- proc.time()
    loglik <- buildCachedObjectiveFun(loglik, MODEL)
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
      control$REPORT <- 0
    }
    do_dream <- function(pars) {
      names(pars) <- names(parlist) ## TODO: dream should retain names
      thisMod <- update(MODEL, newpars = pars)
      if (!isValidModel(thisMod)) {
        return(-1e8)
      }
      obj <- objFunVal(thisMod, objective = loglik)
      if (!is.null(save)) save(pars, obj, thisMod)
      obj
    }
    ans <- dream::dream(do_dream,
      pars = parlist,
      func.type = "logposterior.density",
      control = control
    )
    environment(ans$call) <- environment()
    bestPars <- coef(ans, method = "sample.ml")
    bestModel <- update(MODEL, newpars = bestPars)
    bestModel$funevals <- ans$fun.evals
    bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
    bestModel$objective <- loglik
    if (vcov) {
      ## estimate covariance matrix from final population
      start <- end(ans$Sequences) / 2 + 1
      bestModel$cov.mat <-
        cov(as.matrix(window(ans$Sequences, start = start)))
    }
    bestModel$fit.call <- match.call()
    bestModel$fit.result <- ans
    return(bestModel)
  }
