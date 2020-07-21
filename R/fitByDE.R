## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##




#' Fit a hydromad model using the DE (Differential Evolution) algorithm.
#'
#' Fit a hydromad model using the DE (Differential Evolution) algorithm.
#'
#'
#' @param MODEL a model specification created by \code{\link{hydromad}}. It
#' should not be fully specified, i.e one or more parameters should be defined
#' by \emph{ranges} of values rather than exact values.
#' @param objective objective function to maximise, given as a
#' \code{function(Q, X, ...)}.  See \code{\link{objFunVal}}.
#' @param control settings for the DE algorithm. See
#' \code{\link[DEoptim]{DEoptim.control}}.
#' @return the best model from those sampled, according to the given
#' \code{objective} function. Also, these extra elements are inserted:
#' \item{fit.result}{ the result from \code{\link[DEoptim]{DEoptim}}.  }
#' \item{objective}{ the \code{objective} function used.  } \item{funevals}{
#' total number of evaluations of the model simulation function.  }
#' \item{timing}{ timing vector as returned by \code{system.time}.  }
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link[DEoptim]{DEoptim}}, \code{\link{objFunVal}}
#' @keywords optimization
#' @examples
#'
#' library("DEoptim")
#'
#' data(Cotter)
#' x <- Cotter[1:1000]
#'
#' ## IHACRES CWI model with power law unit hydrograph
#' modx <- hydromad(x, sma = "cwi", routing = "powuh")
#' modx
#'
#' foo <- fitByDE(modx, control = DEoptim.control(itermax = 5))
#'
#' summary(foo)
#'
#' ## return value from DE:
#' str(foo$fit.result)
#'
#' ## plot objective function value convergence over time
#' xyplot(optimtrace(foo),
#'   type = "b",
#'   xlab = "function evaluations", ylab = "objective fn. value"
#' )
#' @export
fitByDE <-
  function(MODEL,
           objective = hydromad.getOption("objective"),
           control = hydromad.getOption("de.control")) {
    if (!requireNamespace("DEoptim")) stop("package DEoptim is required for fitByDE")
    control <- do.call(DEoptim::DEoptim.control, control)
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
    if (!isTRUE(hydromad.getOption("trace"))) {
      control$trace <- FALSE
    }
    lower <- sapply(parlist, min)
    upper <- sapply(parlist, max)
    bestModel <- MODEL
    bestFunVal <- -Inf
    do_de <- function(pars) {
      names(pars) <- names(parlist)
      thisMod <- update(MODEL, newpars = pars)
      if (!isValidModel(thisMod)) {
        return(1e8)
      }
      thisVal <- objFunVal(thisMod, objective = objective)
      if (isTRUE(thisVal > bestFunVal)) {
        bestModel <<- thisMod
        bestFunVal <<- thisVal
      }
      ## DEoptim does minimisation, so:
      return(-thisVal)
    }
    ans <- DEoptim::DEoptim(do_de,
      lower = lower, upper = upper,
      control = control
    )
    bestModel$funevals <- ans$optim$nfeval
    bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
    bestModel$objective <- objective
    bestModel$fit.call <- match.call()
    bestModel$fit.result <- ans
    return(bestModel)
  }
