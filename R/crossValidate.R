## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Joseph Guillaume <josephguillaume@gmail.com>
##



#' Cross-validation of hydromad model specification
#'
#' Using a split dataset, estimate parameters using optimisation on data
#' subsets and evaluate their performance on all other subsets.
#'
#'
#' @param MODEL an object of class \code{hydromad}.
#' @param periods named list of start and end dates, passed to
#' \code{\link{splitData}}
#' @param name.Model.str Name to give to this model structure to allow a
#' combined analysis
#' @param name.Cal.objfn Name to give to this model identification process
#' (e.g. name of objective function and/or optimisation algorithm), to allow a
#' combined analysis
#' @param name.Catchment Name to give to this catchment to allow a combined
#' analysis
#' @param fitBy function to estimate parameters of \code{MODEL}, e.g.
#' \code{\link{fitByOptim}}, \code{\link{fitBySCE}}
#' @param \dots Arguments passed to \code{fitBy}
#' @param trace Whether to report messages.
#' @param parallel name of method to use for parallelisation ("foreach" or
#' "none"), or list giving settings for parallelisation. See
#' \link{hydromad_parallelisation}.
#' @return A runlist of n*n models evaluated in each of n periods with
#' parameters estimated from each of n periods, of subclass
#' \code{crossvalidation}
#' @section Parallelisation: \code{crossValidate} optionally allows the
#' separate optimisations to be run concurrently with the parallel option
#' \code{method="foreach"}. This is usually only worthwhile for longer running
#' optimisations such as \code{\link{fitBySCE}} rather than relatively fast
#' methods such as \code{\link{fitByOptim}}.
#'
#' The total runtime is limited by the runtime of the slowest of optimisations,
#' e.g. the longest data subset, most complex objective function response
#' surface, or slowest computer on which an optimisation is being run. Some of
#' the workers may therefore be idle (potentially wasting money) even though
#' others are still running.
#'
#' The evaluation of parameters on validation data subsets is also optionally
#' parallelised through the function \code{update.runlist}, by setting
#' \code{hydromad.options(parallel=list(update.runlist="clusterApply"))}. The
#' advantage of this is likely to be minor, unless a large number of cross
#' validation periods are used, due to the overhead involved and its relative
#' speed compared to the optimisation. Note that this requires parallelisation
#' to be setup on the worker, which is where the evaluation occurs.
#'
#' If the parallelisation backend for \code{foreach} supports it, the
#' cross-validations can be set to occur in the background using
#' \code{parallel=list(method="foreach",async=TRUE)}. In this case, the
#' function returns immediately and the progress and results can be retrieved
#' using functions provided by the parallelisation backend. This can be useful
#' to submit a number of cross-validations for which the results are not
#' immediately needed. As with a single cross-validation, mixing long and short
#' running optimisations can make it difficult to fully utilise available
#' workers.
#'
#' In future, it may also be possible to parallelise each optimisation itself
#' in addition to or instead of parallelising the optimisation of each data
#' period.
#' @author Joseph Guillaume
#' @seealso \code{\link{paretoTimeAnalysis}}
#' @examples
#'
#' data(Cotter)
#' modx <- hydromad(Cotter,
#'   sma = "cwi", routing = "expuh",
#'   tau_s = c(2, 100), v_s = c(0, 1)
#' )
#' periods <- list(
#'   P1 = as.Date(c("1966-05-01", "1976-04-30")),
#'   P2 = as.Date(c("1976-05-1", "1986-04-30"))
#' )
#' ## Estimate parameters using single fitByOptim run
#' ## from single initial parameter set
#' runs <- crossValidate(modx, periods = periods, fitBy = fitByOptim, samples = 1)
#' summary(runs)
#' ## Cross-validation statistics can then be analysed with other methods
#' paretoTimeAnalysis_areModelsDominated(summary(runs))
#' cast(
#'   melt(summary(runs), id.vars = c("calib.period", "sim.period")),
#'   calib.period ~ variable + sim.period
#' )
#' paretoTimeAnalysis(runs)
#' @export
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
      if (trace) message(sprintf("Running crossvalidate in parallel using foreach: %s", foreach::getDoParName()))
      opts <- hydromad.options()
      runs <- foreach::foreach(
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


#' @export
summary.crossvalidation <- function(object, ...) {
  s <- NextMethod(object, ...)
  s$sim.period <- sapply(object, function(x) x$name.sim.period)
  s$calib.period <- sapply(object, function(x) x$name.calib.period)
  s$Model.str <- sapply(object, function(x) x$name.Model.str)
  s$Cal.objfn <- sapply(object, function(x) x$name.Cal.objfn)
  s$Catchment <- sapply(object, function(x) x$name.Catchment)
  s
}

#' @import utils
utils::globalVariables(c("%dopar%"))
