#' Multi-objective optimisation using NSGAII
#'
#' Estimate multi-objective Pareto front using NSGAII
#'
#'
#' @param MODEL a model specification created by \code{\link{hydromad}}. It
#' should not be fully specified, i.e one or more parameters should be defined
#' by \emph{ranges} of values rather than exact values.
#' @param objective objective functions to maximise, as a list with elements as
#' \code{function(Q, X, ...)}.  See \code{\link{objFunVal}}.
#' @param control arguments for nsga2 function. See \code{nsga2}.
#' @return \code{\link{runlist}} of models on Pareto front
#' @author Joseph Guillaume
#' @keywords optimization
#' @examples
#'
#' data(Cotter)
#' x <- Cotter[1:1000]
#'
#' ## IHACRES CWI model with exponential unit hydrograph
#' ## an unfitted model, with ranges of possible parameter values
#' modx <- hydromad(x,
#'   sma = "cwi", routing = "expuh",
#'   tau_s = c(2, 100), v_s = c(0, 1)
#' )
#' ## Multi-objective optimisation
#' front <- paretoObjectivesNsga2(modx, objective = list(hmadstat("r.sq.log"), hmadstat("r.squared")))
#' ## Pairwise plot of parameters on Pareto front
#' splom(coef(front))
#' ## Calculate objectives
#' stats <- t(sapply(front, objFunVal, objective = list(hmadstat("r.sq.log"), hmadstat("r.squared"))))
#' splom(stats)
#' @export
paretoObjectivesNsga2 <- function(MODEL, objective = hydromad.getOption("objective"),
                                  control = hydromad.getOption("nsga2.control")) {
  if (!requireNamespace("mco")) stop("package mco is required for paretoObjectivesNsga2")
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
  lower <- sapply(parlist, min)
  upper <- sapply(parlist, max)
  do_nsga2 <- function(pars) {
    names(pars) <- names(parlist)
    thisMod <- update(MODEL, newpars = pars)
    if (!isValidModel(thisMod)) {
      return(1e+08)
    }
    thisVal <- objFunVal(thisMod, objective = objective)
    return(-unlist(thisVal))
  }
  args <- modifyList(control, list(
    fn = do_nsga2,
    idim = length(parlist),
    odim = length(objective),
    lower.bounds = lower, upper.bounds = upper
  ))
  ans <- do.call(mco::nsga2, args)
  colnames(ans$par) <- names(parlist)
  MODEL$funevals <- NA ## TODO
  MODEL$timing <- signif(proc.time() - start_time, 4)[1:3]
  MODEL$objective <- objective
  MODEL$fit.call <- match.call()
  MODEL$fit.result <- ans
  if (isTRUE(hydromad.getOption("trace"))) cat("Creating runlist (running models on pareto front)\n")
  ## TODO: potential for parallelisation
  front <- as.runlist(apply(ans$par, 1, function(p) update(MODEL, newpars = p)))
  return(front)
}
