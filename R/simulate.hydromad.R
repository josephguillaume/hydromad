## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##




#' Simulate hydromad models by parameter sampling.
#'
#' Run many simulations by sampling within parameter ranges.
#'
#' none yet.
#'
#' @importFrom stats simulate
#'
#' @param object a \code{hydromad} object (produced by the
#' \code{\link{hydromad}()} function) that is not fully specified (i.e. so some
#' parameter values are given as ranges).
#' @param nsim number of parameter samples to run.
#' @param seed optional random seed, for repeatability.
#' @param \dots further arguments to \code{FUN}.
#' @param sampletype sampling method; see \code{\link{parameterSets}}.
#' @param FUN optional function to apply to each simulated model. Typical
#' examples would be \code{\link{objFunVal}},
#' \code{\link[=summary.hydromad]{summary}} or
#' \code{\link[=predict.hydromad]{predict}}.
#' @param objective an objective function (or statistic function); this is just
#' an argument to be passed on to \code{FUN}, which in this case defaults to
#' \code{\link{objFunVal}} to calculate the statistic value from each model. It
#' is treated as a special argument because it is cached before the simulation
#' run for efficiency.
#' @param bind to bind the result from \code{FUN} as one or more columns onto
#' the matrix of parameter values, and return as a data frame.
#' @return a list of results, where each element is named by its parameter set.
#' The result also has an attribute \code{"psets"} which gives the parameter
#' values used in each simulation (as a data frame).
#'
#' If \code{bind = TRUE}, a data frame.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{parameterSets}}, \code{\link{fitBySampling}}
#' @keywords iteration
#' @examples
#'
#' data(Canning)
#' mod0 <- hydromad(Canning[1:500, ], sma = "cwi")
#' sim0 <- simulate(mod0, nsim = 5, sampletype = "latin")
#' coef(sim0)
#' summary(sim0)
#'
#' ## plot the objective function surface over two parameters
#' mod1 <- update(mod0, routing = "armax", rfit = list("ls", order = c(2, 1)))
#' sim1 <- simulate(mod1, 144,
#'   sampletype = "all", FUN = objFunVal,
#'   objective = ~ nseStat(Q, X, trans = sqrt)
#' )
#' levelplot(result ~ tw + f, sim1,
#'   cex = 2,
#'   panel = panel.levelplot.points,
#'   main = "R Squared (of sq.rt. data) over parameter space"
#' ) +
#'   latticeExtra::layer(panel.2dsmoother(...), under = TRUE)
#'
#' ## dotty plots (list any number of parameters in formula)
#' xyplot(result ~ tw + f, sim1, outer = TRUE)
#' @export
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
