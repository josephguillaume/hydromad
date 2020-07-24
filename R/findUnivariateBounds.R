#' Find univariate feasible bounds
#'
#' Find bounds on parameters within which a threshold of performance is
#' satisfied by solving \code{objective(single.par)=thres} using univariate
#' line searches for the lower and upper bounds of each parameter.
#'
#' This is intended as a heuristic method for identifying feasible bounds. It
#' may not work as desired in some cases.
#'
#' Errors may be produced by \code{\link{uniroot}} of the type \code{ f()
#' values at end points not of opposite sign}. This indicates that no parameter
#' value could be found that has an objective value equal to \code{thres}. In
#' that case, the relevant bound is kept unchanged as specified in \code{modx}.
#'
#' Errors of the type \code{lower < upper is not fulfilled} occur if the
#' parameters of \code{fitx} are at the boundaries specified in \code{modx}.
#' This error can usually be ignored, though it suggests that the \code{fitx}
#' parameter estimation could perhaps be improved.
#'
#' @importFrom stats uniroot
#'
#' @param modx a model specification created by \code{\link{hydromad}}. It
#' should not be fully specified, i.e one or more parameters should be defined
#' by \emph{ranges} of values rather than exact values.
#' @param fitx a best guess model with fully specified parameters. This will be
#' used as one end point for a line search
#' @param thres Minimum threshold on objective, defining a feasible region.
#' Should be less than \code{objFunVal(fitx,objective=objective)}
#' @param objective objective function used to determine feasible region, given
#' as a \code{function(Q, X, ...)}.  See \code{\link{objFunVal}}.
#' @return A named list of bounds for each free parameter in \code{modx}.
#' @author Joseph Guillaume
#' @seealso \code{\link{runRSM}} which requires narrowed bounds
#' @keywords models
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
#'
#' ## Best fit, used as initial feasible solution
#' fitx <- fitByOptim(modx)
#'
#' ## Identify bounds
#' thres <- objFunVal(fitx) - 0.05 ## Look for models within 0.05 of best fit
#' bounds <- findUnivariateBounds(modx, fitx, thres)
#'
#' ## Create model object with the new bounds
#' modx2 <- update(modx, newpars = bounds)
#'
#' bounds
#' modx2
#' @export
findUnivariateBounds <- function(modx, fitx, thres, objective = hydromad.getOption("objective")) {
  pars <- getFreeParsRanges(modx)
  for (p in names(pars)) {
    best <- coef(fitx)[[p]]
    ## Lower side
    try({
      lo <- uniroot(function(val) {
        names(val) <- p
        objFunVal(update(fitx, newpars = val), objective = objective) - thres
      }, interval = c(min(pars[[p]]), best))
      pars[[p]][1] <- lo$root
    })
    ## upper side
    try({
      hi <- uniroot(function(val) {
        names(val) <- p
        objFunVal(update(fitx, newpars = val), objective = objective) - thres
      }, interval = c(best, max(pars[[p]])))
      pars[[p]][2] <- hi$root
    })
  }
  return(pars)
}
