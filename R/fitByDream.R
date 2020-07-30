## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##



#' Fit a hydromad model using the DREAM (DiffeRential Evolution Adaptive
#' Metropolis) algorithm.
#'
#' Fit a hydromad model using the DREAM (DiffeRential Evolution Adaptive
#' Metropolis) algorithm. This is a Markov Chain Monte Carlo algorithm which
#' gives estimates of the joint probability distribution of parameters
#' according to a likelihood function. The fitting function returns the maximum
#' likelihood model, but the full MCMC results are also available as component
#' \code{$fit.result}. The result can also be used to define a
#' \link[=defineFeasibleSet]{feasible parameter set}.
#'
#' @importFrom stats cov
#'
#' @param MODEL a model specification created by \code{\link{hydromad}}. It
#' should not be fully specified, i.e one or more parameters should be defined
#' by \emph{ranges} of values rather than exact values.
#' @param loglik log-likelihood function (log of the posterior probability
#' density), given as a \code{function(Q, X, ...)}.  See
#' \code{\link{objFunVal}}.
#' @param control settings for the DREAM algorithm. See
#' \code{\link[dream]{dream}}.
#' @param vcov if \code{vcov = TRUE}, the parameter variance-covariance matrix
#' will be estimated from the last half of the sequences.  It can be extract
#' using \code{\link{vcov}}.
#' @param save Optional \code{function(pars,objective value,model)} that will
#' be called for every model evaluation, for example to save every model run.
#' @return the best model from those sampled, according to the given
#' \code{loglik} function. Also, these extra elements are inserted:
#' \item{fit.result}{ the result from \code{\link[dream]{dream}}.  }
#' \item{objective}{ the \code{loglik} function used.  } \item{funevals}{ total
#' number of evaluations of the model simulation function.  } \item{timing}{
#' timing vector as returned by \code{system.time}.  }
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link[dream]{dream}}, \code{\link{objFunVal}}
#' @keywords optimization
#' @examples
#'
#'
#' if (requireNamespace("dream", quietly = TRUE)) {
#'   data(Cotter)
#'   x <- Cotter[1:1000]
#'
#'   ## IHACRES CWI model with power law unit hydrograph
#'   modx <- hydromad(x, sma = "cwi", routing = "powuh")
#'   modx
#'
#'   ## a very short run! just to demonstrate methods
#'   foo <- fitByDream(modx, control = list(ndraw = 500))
#'
#'   summary(foo)
#'
#'   ## parameter correlation matrix with symbols
#'   symnum(cov2cor(vcov(foo)))
#'
#'   ## return value from dream:
#'   str(foo$fit.result)
#'
#'   ## plot log-likelihood value convergence over time
#'   xyplot(window(optimtrace(foo, raw = TRUE), start = 50),
#'     superpose = TRUE, auto.key = FALSE,
#'     xlab = "function evaluations", ylab = "neg. log likelihood"
#'   )
#'
#'   ## calculate corresponding objective function values over time.
#'   xyplot(optimtrace(foo, objective = ~ -hmadstat("r.squared")(Q, X)),
#'     xlab = "function evaluations", ylab = "negative R Squared"
#'   )
#'
#'   ## MCMC diagnostics and more are available:
#'   methods(class = "dream")
#' }
#' @export
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
