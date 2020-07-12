## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Joseph Guillaume <josephguillaume@gmail.com>
##

crossValidate <- function(MODEL, periods,
                          name.Model.str = paste(MODEL$sma, MODEL$routing),
                          name.Cal.objfn = "unknown",
                          name.Catchment = as.character(MODEL$call$DATA),
                          fitBy, ...,
                          trace = isTRUE(hydromad.getOption("trace")),
                          parallel = hydromad.getOption("parallel")[["crossValidate"]]) {
  cv.set <- splitData(MODEL, periods = periods)

  # Sets default settings for parallelisation if missing
  parallel <- hydromad.parallel(parallel)

  switch(parallel$method,
    "foreach" = {
      if (trace) message(sprintf("Running crossvalidate in parallel using foreach: %s", getDoParName()))
      opts <- hydromad.options()
      runs <- foreach(
        n = names(cv.set),
        .packages = parallel$packages,
        .inorder = FALSE,
        .export = parallel$export,
        .final = function(runs) {
          runs <- do.call(c, runs)
          class(runs) <- unique(c("crossvalidation", class(runs), "runlist", "list"))
          runs
        },
        .options.redis = list(async = parallel$async)
      ) %dopar% {
        hydromad.options(opts)
        if (trace) cat("\nFitting period: ", n, "\n")
        fitx <- fitBy(cv.set[[n]], ...)
        new.runs <- update(cv.set, newpars = coef(fitx))
        names(new.runs) <- sprintf("%s_cal%s", names(cv.set), n)
        ## Preserve fit attributes
        new.runs[[sprintf("%s_cal%s", n, n)]] <- fitx
        for (m in 1:length(new.runs)) {
          new.runs[[m]]$name.Model.str <- name.Model.str
          new.runs[[m]]$name.Cal.objfn <- name.Cal.objfn
          new.runs[[m]]$name.calib.period <- n
          new.runs[[m]]$name.sim.period <- names(cv.set)[m]
          new.runs[[m]]$name.Catchment <- name.Catchment
        }
        new.runs
      }
    },
    {
      if (trace) message("Running crossvalidate sequentially")
      runs <- runlist()
      for (n in names(cv.set)) {
        if (trace) cat("\nFitting period: ", n, "\n")
        fitx <- fitBy(cv.set[[n]], ...)
        new.runs <- update(cv.set, newpars = coef(fitx))
        names(new.runs) <- sprintf("%s_cal%s", names(cv.set), n)
        ## Preserve fit attributes
        new.runs[[sprintf("%s_cal%s", n, n)]] <- fitx
        for (m in 1:length(new.runs)) {
          new.runs[[m]]$name.Model.str <- name.Model.str
          new.runs[[m]]$name.Cal.objfn <- name.Cal.objfn
          new.runs[[m]]$name.calib.period <- n
          new.runs[[m]]$name.sim.period <- names(cv.set)[m]
          new.runs[[m]]$name.Catchment <- name.Catchment
        }
        runs <- c(runs, new.runs)
      }
      class(runs) <- unique(c("crossvalidation", class(runs), "runlist", "list"))
    }
  ) ## switch parallel
  return(runs)
}

summary.crossvalidation <- function(object, ...) {
  s <- NextMethod(object, ...)
  s$sim.period <- sapply(object, function(x) x$name.sim.period)
  s$calib.period <- sapply(object, function(x) x$name.calib.period)
  s$Model.str <- sapply(object, function(x) x$name.Model.str)
  s$Cal.objfn <- sapply(object, function(x) x$name.Cal.objfn)
  s$Catchment <- sapply(object, function(x) x$name.Catchment)
  s
}
