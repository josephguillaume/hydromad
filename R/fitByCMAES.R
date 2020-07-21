#' Fit a hydromad model using CMA-ES (Covariance matrix adapting evolutionary
#' strategy) from cmaes package
#'
#'
#' Placeholder for description
#'
#' @param MODEL a model specification created by \code{\link{hydromad}}. It
#' should not be fully specified, i.e one or more parameters should be defined
#' by \emph{ranges} of values rather than exact values.
#' @param objective objective function to maximise, given as a
#' \code{function(Q, X, ...)}.  See \code{\link{objFunVal}}.
#' @param control settings for the CMA-ES algorithm. See \code{cma_es}
#' @param vcov Ignored
#' @return the best model from those sampled, according to the given
#' \code{objective} function. Also, these extra elements are inserted:
#' \item{fit.result}{ the result from \code{\link{SCEoptim}}.  }
#' \item{objective}{ the \code{objective} function used.  } \item{funevals}{
#' total number of evaluations of the model simulation function.  }
#' \item{timing}{ timing vector as returned by \code{system.time}.  }
#' @author Joseph Guillaume
#' @seealso \code{cma_es}
#' @references Placeholder
#' @keywords optimisation
#' @export
fitByCMAES <- function(MODEL, objective = hydromad.getOption("objective"),
                       control = hydromad.getOption("cmaes.control"), vcov = FALSE) {
  if (!requireNamespace("cmaes")) stop("package cmaes is required for fitByCMAES")
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
  ## TODO: if (isTRUE(hydromad.getOption("trace")))
  lower <- sapply(parlist, min)
  upper <- sapply(parlist, max)
  initpars <- sapply(parlist, mean)
  bestModel <- MODEL
  bestFunVal <- -Inf
  do_cmaes <- function(pars) {
    thisMod <- update(MODEL, newpars = pars)
    if (!isValidModel(thisMod)) {
      return(NA)
    }
    thisVal <- objFunVal(thisMod, objective = objective)
    if (isTRUE(thisVal > bestFunVal)) {
      bestModel <<- thisMod
      bestFunVal <<- thisVal
    }
    return(-thisVal)
  }
  ans <- cmaes::cma_es(initpars, do_cmaes,
    lower = lower, upper = upper,
    control = control
  )
  if (ans$convergence != 0) {
    if (!isTRUE(hydromad.getOption("quiet"))) {
      warning(ans$message)
    }
    bestModel$msg <- ans$message
  }
  bestModel$funevals <- ans$counts[1]
  bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
  bestModel$objective <- objective
  if (vcov) {
    warning("vcov not yet implemented")
  }
  bestModel$fit.call <- match.call()
  bestModel$fit.result <- ans
  return(bestModel)
}
