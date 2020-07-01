## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


fitByOptim <-
  function(MODEL,
           objective = hydromad.getOption("objective"),
           method = hydromad.getOption("optim.method"), control = list(),
           samples = hydromad.getOption("fit.samples"),
           sampletype = c("latin.hypercube", "random", "all.combinations"),
           initpars = NULL, multistart = FALSE,
           vcov = FALSE, hessian = vcov) {
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
    ## if only one parameter, use specialised function
    if (length(parlist) == 1) {
      warning("only one parameter; switching to fitByOptim1")
      return(fitByOptim1(MODEL, objective = objective))
    }
    if (vcov && (method == "PORT")) {
      warning("vcov does not work with method = 'PORT'")
      vcov <- FALSE
    }
    if (multistart) {
      ## multiple optimisation runs with different starting points.
      ## generate parameter sets
      psets <- parameterSets(parlist, samples = samples, method = sampletype)
      bestModel <- MODEL
      bestFunVal <- -Inf
      objseq <- numeric()
      funevals <- 0
      for (i in seq(NROW(psets))) {
        thisPars <- as.list(psets[i, , drop = FALSE])
        thisMod <-
          try(fitByOptim(MODEL,
            objective = objective,
            method = method, control = control,
            hessian = hessian,
            multistart = FALSE, initpars = thisPars
          ),
          silent = TRUE
          )
        if (!isValidModel(thisMod)) {
          next
        }
        funevals <- funevals + thisMod$funevals
        objseq <- c(objseq, thisMod$fit.result$objseq)
        thisVal <- objFunVal(thisMod, objective = objective)
        if (thisVal > bestFunVal) {
          bestModel <- thisMod
          bestFunVal <- thisVal
        }
      }
      bestModel$funevals <- funevals
      bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
      bestModel$fit.call <- match.call()
      bestModel$fit.result$objseq <- objseq
      return(bestModel)
    } else {
      ## single optimisation run.
      lower <- sapply(parlist, min)
      upper <- sapply(parlist, max)
      #        lower[is.na(lower)] <- -Inf
      #        upper[is.na(upper)] <- Inf
      pre.funevals <- 0
      objseq <- numeric()
      preMODEL <- MODEL
      if (!is.null(initpars)) {
        ## initial parameter values were specified.
        initpars <- unlist(initpars)[names(parlist)]
      } else if (samples >= 1) {
        ## do sampling to find a good starting point
        preMODEL <- fitBySampling(MODEL,
          objective = objective,
          samples = samples, sampletype = sampletype
        )
        if (!isValidModel(preMODEL)) {
          return(preMODEL)
        }
        pre.funevals <- preMODEL$funevals
        objseq <- preMODEL$fit.result$objseq
        initpars <- coef(preMODEL)[names(parlist)]
      } else {
        initpars <- sapply(parlist, mean)
      }
      if (!identical(names(initpars), names(parlist))) stop("Names of initpars do not match parlist. If using rfit, check that order matches the parameters specified to hydromad")
      ## now optimise
      control0 <- hydromad.getOption("optim.control")
      ## choose parscale automatically (nlminb has 'scale' = 1/parscale)
      parscale <- control0$parscale
      if (is.null(parscale)) {
        parscale <- ifelse(initpars == 0, 1, abs(initpars))
      }
      control0$parscale <- parscale
      if (method == "PORT") {
        control0 <- hydromad.getOption("nlminb.control")
      }
      control <- modifyList(control0, control)
      if (isTRUE(hydromad.getOption("quiet"))) {
        control$trace <- 0
      }
      bestModel <- MODEL
      bestFunVal <- -Inf
      i <- length(objseq)
      objseq <- c(objseq, rep(NA_real_, 100))
      do_optim <- function(pars) {
        ## TODO: handle boundaries better
        ## just set params to bounded values?
        k <- which(pars < lower)
        if (any(k)) {
          return(1e12 + (sum(lower[k] - pars[k])) * 1e6)
        }
        k <- which(pars > upper)
        if (any(k)) {
          return(1e12 + (sum(pars[k] - upper[k])) * 1e6)
        }
        # if (any(pars < lower)) return(NA)
        # if (any(pars > upper)) return(NA)

        i <<- i + 1

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
        ## optim does minimisation, so:
        return(-thisVal)
      }
      if (!isTRUE(hydromad.getOption("catch.errors.optim"))) {
        try <- force
      } ## i.e. skip the try()
      lowerb <- -Inf
      upperb <- Inf
      if (method == "L-BFGS-B") {
        lowerb <- lower
        upperb <- upper
      }
      if (method == "PORT") {
        ans <- try(nlminb(initpars, do_optim,
          lower = lower, upper = upper,
          control = control, scale = 1 / parscale
        ))
      } else {
        ans <- try(optim(initpars, do_optim,
          method = method,
          lower = lowerb, upper = upperb,
          control = control, hessian = hessian
        ))
      }
      if (inherits(ans, "try-error")) {
        bestModel$msg <- ans
        return(bestModel)
      }
      if (ans$convergence != 0) {
        msg <- if (method == "PORT") {
          ans$message
        } else if (ans$convergence == 1) {
          "optim() reached maximum iterations"
        } else {
          paste(
            "optim() returned convergence code",
            ans$convergence
          )
        }
        if (!isTRUE(hydromad.getOption("quiet"))) {
          warning(msg)
        }
        bestModel$msg <- msg
      }
      bestModel$funevals <- i
      bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
      bestModel$objective <- objective
      if (vcov) {
        ## approximate covariance matrix from inverse of hessian (often poor!)
        bestModel$cov.mat <- solve(ans$hessian)
      }
      bestModel$fit.call <- match.call()
      ans$objseq <- objseq[1:i]
      bestModel$fit.result <- ans
      return(bestModel)
    }
  }
