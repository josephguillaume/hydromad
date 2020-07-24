## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##




#' Fit a hydromad model using general optimisation algorithms.
#'
#' Fits a hydromad model using R's \code{optim} or \code{nlminb} functions. Has
#' multi-start and pre-sampling options.
#'
#' See \code{\link{optim}} for a brief description of the available algorithms,
#' including references to literature.
#'
#' \code{\link{fitByOptim1}} handles a single free parameter using
#' \code{\link{optimize}}; it is called by \code{fitByOptim}, with a warning,
#' in that case.
#'
#' @importFrom stats optim
#'
#' @aliases fitByOptim fitByOptim1
#' @param MODEL a model specification created by \code{\link{hydromad}}. It
#' should not be fully specified, i.e one or more parameters should be defined
#' by \emph{ranges} of values rather than exact values.
#' @param objective objective function to maximise, given as a
#' \code{function(Q, X, ...)}.  See \code{\link{objFunVal}}.
#' @param method,control optimisation algorithm and settings. See
#' \code{\link{optim}}. The default is \code{"Nelder-Mead"} (a simplex
#' algorithm).  An additional method \code{"PORT"} can be given, which calls
#' \code{\link{nlminb}}.
#' @param samples number of parameter sets to test.
#' @param sampletype sampling scheme -- see \code{\link{parameterSets}}.
#' @param initpars initial parameter set.
#' @param multistart if this is \code{TRUE}, then each of the initial parameter
#' samples (based on \code{samples} and \code{sampletype}) are used as a
#' starting value for a separate optimisation run, and the best result is kept.
#' @param vcov,hessian if \code{vcov = TRUE}, the parameter variance-covariance
#' matrix will be estimated as the inverse of the hessian matrix returned by
#' \code{optim()}. This may be a very poor estimate! It can be extract using
#' \code{\link{vcov}}.  Does not work with \code{method = "PORT"}.
# @param tol convergence tolerance, see \code{\link{optimize}}.
#' @return the best model from those sampled, according to the given
#' \code{objective} function. Also, these extra elements are inserted:
#' \item{fit.result}{ the result from \code{\link{optim}}.  } \item{objective}{
#' the \code{objective} function used.  } \item{funevals}{ total number of
#' evaluations of the model simulation function.  } \item{timing}{ timing
#' vector as returned by \code{system.time}.  }
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{fitBySampling}}, \code{\link{optim}},
#' \code{\link{objFunVal}}
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
#' set.seed(0)
#' foo <- fitByOptim(modx, samples = 100)
#'
#' summary(foo)
#'
#' ## return value from optim (with 'objseq' added):
#' str(foo$fit.result)
#'
#' ## plot objective function value convergence over time
#' xyplot(optimtrace(foo),
#'   type = "b",
#'   xlab = "function evaluations", ylab = "objective fn. value"
#' )
#'
#' ## repeat optimisation with single random starting points
#' fooreps <-
#'   replicate(4,
#'     fitByOptim(modx,
#'       samples = 1, sampletype = "random",
#'       objective = hmadstat("r.squared"),
#'       method = "Nelder-Mead"
#'     ),
#'     simplify = FALSE
#'   )
#' names(fooreps) <- paste("rep.", seq_along(fooreps))
#' ## extract and plot the optimisation traces
#' traces <- lapply(fooreps, optimtrace)
#' tracesraw <- lapply(fooreps, optimtrace, raw = TRUE)
#' xyplot(do.call("merge", traces),
#'   superpose = TRUE,
#'   sub = 'method = "Nelder-Mead"',
#'   xlab = "Fn. evaluations", ylab = "Objective value",
#'   auto.key = list(corner = c(1, 0))
#' ) +
#'   xyplot(do.call("merge", tracesraw),
#'     superpose = TRUE,
#'     type = "p", cex = 0.5
#'   )
#'
#' ## if you try it again with method = "PORT" you will find that all
#' ## replicates converge to the optimum regardless of starting point.
#' @export
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
