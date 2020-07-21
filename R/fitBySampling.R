## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##




#' Fit a hydromad model by sampling the parameter space.
#'
#' Fit a hydromad model by sampling the parameter space.  Returns best result
#' from sampling in parameter ranges using random, latin hypercube sampling, or
#' a uniform grid (all combinations). The function also retains the parameter
#' sets and objective function values, which can be used to define a
#' \link[=defineFeasibleSet]{feasible parameter set}
#'
#' See \code{\link{parameterSets}}.
#'
#' @param MODEL a model specification created by \code{\link{hydromad}}. It
#' should not be fully specified, i.e one or more parameters should be defined
#' by \emph{ranges} of values rather than exact values.
#' @param objective objective function to maximise, given as a
#' \code{function(Q, X, ...)}.  See \code{\link{objFunVal}}.
#' @param samples number of parameter sets to test.
#' @param sampletype sampling scheme -- see \code{\link{parameterSets}}.
#' @return the best model from those sampled, according to the given
#' \code{objective} function.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{fitByOptim}}, \code{\link{parameterSets}},
#' \code{\link{objFunVal}}
#' @keywords optimization
#' @examples
#'
#'
#' data(Cotter)
#' x <- Cotter[1:1000]
#'
#' ## IHACRES CWI model with armax unit hydrograph fitted by least squares
#' modx <- hydromad(x, sma = "cwi", routing = "armax", rfit = "ls")
#' modx
#'
#' foo <- fitBySampling(modx)
#'
#' summary(foo)
#'
#' ## plot objective function value improvement over time
#' xyplot(optimtrace(foo),
#'   type = "b",
#'   xlab = "function evaluations", ylab = "objective fn. value"
#' )
#' @export
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
