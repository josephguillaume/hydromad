## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


simulate.hydromad <-
  function(object, nsim, seed, ...,
           sampletype = c("latin.hypercube", "random", "all.combinations"),
           FUN = NULL, objective = NULL, bind = !is.null(objective)) {
    sampletype <- match.arg(sampletype)
    MODEL <- object
    samples <- nsim
    if (!missing(seed)) {
      set.seed(seed)
    }
    if (is.character(FUN)) {
      FUN <- get(FUN, mode = "function")
    }
    if (!is.null(objective)) {
      ## catch the 'objective' argument, cache it, and pass on to FUN
      objective <- buildCachedObjectiveFun(objective, MODEL)
      if (is.null(FUN)) FUN <- objFunVal
      origFUN <- FUN
      FUN <- function(...) origFUN(..., objective = objective)
    }
    ## this is mostly a copy of fitBySampling
    parlist <- as.list(coef(MODEL, warn = FALSE))
    ## remove any missing parameters
    isok <- sapply(parlist, function(x) !any(is.na(x)))
    parlist <- parlist[isok]
    ## check which parameters are uniquely specified
    isfixed <- (sapply(parlist, length) == 1)
    if (all(isfixed)) {
      warning("all parameters are fixed, so can not fit")
      return(list(if (is.null(FUN)) MODEL else FUN(MODEL, ...)))
    }
    ## generate parameter sets
    psets <- parameterSets(parlist, samples = samples, method = sampletype)
    result <- list()
    length(result) <- NROW(psets)
    for (i in seq_len(NROW(psets))) {
      thisPars <- as.list(psets[i, , drop = FALSE])
      run_name <- paste(names(thisPars), signif(unlist(thisPars), 3),
        sep = "=", collapse = ", "
      )
      if (isTRUE(hydromad.getOption("trace"))) {
        message(run_name)
      }
      thisMod <- update(MODEL, newpars = thisPars)
      ## store model or derived result
      result[[i]] <-
        if (is.null(FUN)) thisMod else FUN(thisMod, ...)
      names(result)[i] <- run_name
    }
    if (bind) {
      if (is.null(FUN)) stop("bind requires the 'FUN' argument")
      names(result) <- NULL
      result <- lapply(result, unlist)
      lens <- range(unlist(lapply(result, length)))
      if (max(lens) != min(lens)) {
        result <- lapply(result, function(x) {
          length(x) <- max(lens)
          x
        })
      }
      result <- do.call(rbind, result)
      return(cbind(psets, result))
    }
    if (is.null(FUN)) {
      result <- as.runlist(result)
    }
    attr(result, "psets") <- psets
    result
  }
