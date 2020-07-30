## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##




#' Fit a hydromad model using the SCE (Shuffled Complex Evolution) algorithm.
#'
#' Fit a hydromad model using the SCE (Shuffled Complex Evolution) algorithm.
#'
#'
#' @param MODEL a model specification created by \code{\link{hydromad}}. It
#' should not be fully specified, i.e one or more parameters should be defined
#' by \emph{ranges} of values rather than exact values.
#' @param objective objective function to maximise, given as a
#' \code{function(Q, X, ...)}.  See \code{\link{objFunVal}}.
#' @param control settings for the SCE algorithm. See \code{\link{SCEoptim}}.
#' Note that unlike SCEoptim, the objective function is maximised
#' (\code{fnscale=-1}) by default.
#' @param vcov if \code{vcov = TRUE}, the parameter variance-covariance matrix
#' will be estimated from the final population.  It can be extract using
#' \code{\link{vcov}}.
#' @return the best model from those sampled, according to the given
#' \code{objective} function. Also, these extra elements are inserted:
#' \item{fit.result}{ the result from \code{\link{SCEoptim}}.  }
#' \item{objective}{ the \code{objective} function used.  } \item{funevals}{
#' total number of evaluations of the model simulation function.  }
#' \item{timing}{ timing vector as returned by \code{system.time}.  }
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{SCEoptim}}, \code{\link{objFunVal}}
#' @keywords optimization
#' @examples
#'
#' data(Cotter)
#' x <- Cotter[1:1000]
#'
#' ## IHACRES CWI model with power law unit hydrograph
#' modx <- hydromad(x, sma = "cwi", routing = "powuh")
#' modx
#'
#' ## run with cut-down settings (for a speedy example only!)
#' foo <- fitBySCE(modx, control = list(maxit = 5, ncomplex = 2))
#'
#' summary(foo)
#'
#' ## return value from SCE:
#' str(foo$fit.result)
#'
#' ## plot objective function value convergence over time
#' xyplot(optimtrace(foo, raw = TRUE),
#'   screens = 1, type = "p",
#'   jitter.x = TRUE, ylim = c(0.7, NA), xlim = c(0, NA),
#'   xlab = "function evaluations", ylab = "objective fn. value"
#' ) +
#'   latticeExtra::layer(panel.average(..., horiz = FALSE, fun = max, lwd = 2))
#' @export
fitBySCE <-
  function(MODEL,
           objective = hydromad.getOption("objective"),
           control = hydromad.getOption("sce.control"),
           vcov = FALSE) {
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
    ## Maximise by default
    control <- modifyList(list(fnscale = -1), control)
    ## remove any fixed parameters
    parlist <- parlist[!isfixed]
    if (isTRUE(hydromad.getOption("trace"))) {
      control$trace <- 1
    }
    lower <- sapply(parlist, min)
    upper <- sapply(parlist, max)
    initpars <- sapply(parlist, mean) ## TODO: allow sampling?
    bestModel <- MODEL
    bestFunVal <- Inf * control$fnscale
    do_sce <- function(pars) {
      thisMod <- update(MODEL, newpars = pars)
      if (!isValidModel(thisMod)) {
        return(NA)
      }
      thisVal <- objFunVal(thisMod, objective = objective)
      if (isTRUE(thisVal * control$fnscale < bestFunVal * control$fnscale)) {
        bestModel <<- thisMod
        bestFunVal <<- thisVal
      }
      ## We use fnscale, so SCEoptim deals with maximisation vs minimisation
      return(thisVal)
    }
    ans <- SCEoptim(do_sce, initpars,
      lower = lower, upper = upper,
      control = control
    )
    if (ans$convergence != 0) {
      if (!isTRUE(hydromad.getOption("quiet"))) {
        warning(ans$message)
      }
      bestModel$msg <- ans$message
    }
    bestModel$funevals <- ans$counts
    bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
    bestModel$objective <- objective
    if (vcov) {
      ## estimate covariance matrix from final population
      ## TODO
      # bestModel$cov.mat <-
      # ans$POP.ALL
      warning("vcov not yet implemented")
    }
    bestModel$fit.call <- match.call()
    bestModel$fit.result <- ans
    return(bestModel)
  }
