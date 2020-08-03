#' Multi-objective optimisation by varying weights
#'
#' Estimate multi-objective Pareto front using multiple weighted single
#' objective optimisations
#'
#' @importFrom parallel clusterApply parApply clusterExport
#'
#' @param MODEL a model specification created by \code{\link{hydromad}}. It
#' should not be fully specified, i.e one or more parameters should be defined
#' by \emph{ranges} of values rather than exact values.
#' @param objective objective functions to maximise, as a list with elements as
#' \code{function(Q, X, ...)}.  See \code{\link{objFunVal}}.
#' @param weights matrix of weights to use, with columns in the same order as
#' \code{objective}
#' @param fitBy function to estimate parameters of \code{MODEL}, e.g.
#' \code{\link{fitByOptim}}, \code{\link{fitBySCE}}
#' @param \dots Arguments passed to \code{fitBy}
#' @return \code{\link{runlist}} of models on Pareto front, corresponding to
#' the \code{weights} specified for each \code{objective}
#' @author Joseph Guillaume
#' @seealso \code{\link{paretoObjectivesNsga2}}
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
#'
#' ## Uncomment to parallelise the fitBy runs
#' # library(parallel)
#' # cl <- makeCluster(getOption("cl.cores", 2))
#' # clusterEvalQ(cl,library(hydromad))
#' # hydromad.options(parallel="clusterApply")
#'
#'
#' ## Optimisation of multiple weights for r.sq.log and r.squared
#' weights <- cbind(c(0, 0.33, 0.5, 0.67, 1), 1 - c(0, 0.33, 0.5, 0.67, 1))
#'
#' ## Estimate parameters using single fitByOptim run
#' ## from single initial parameter set
#' front <- paretoObjectivesVaryWeights(modx,
#'   objective = list(hmadstat("r.sq.log"), hmadstat("r.squared")),
#'   weights = weights, fitBy = fitByOptim, samples = 1
#' )
#'
#' summary(front)
#' ## Plot objectives
#' stats <- t(sapply(front, objFunVal, objective = list(hmadstat("r.sq.log"), hmadstat("r.squared"))))
#' plot(stats)
#' @export
paretoObjectivesVaryWeights <- function(MODEL, objective = hydromad.getOption("objective"), weights, fitBy, ...) {
  objective <- buildCachedObjectiveFun(objective, MODEL)
  switch(hydromad.getOption("parallel")[["paretoObjectivesVaryWeights"]],
    "clusterApply" = {
      if (!requireNamespace("parallel")) stop('package parallel is required for paretoObjectivesVaryWeights if hydromad.getOption("parallel")[["paretoObjectivesVaryWeights"]]=="clusterApply"')
      clusterExport(cl, c("objective", "fitBy", "MODEL"), envir = environment())
      parallel::clusterExport(cl, c("objective", "fitBy", "MODEL"), envir = environment())
      front <- parallel::parApply(
        cl, weights, 1,
        ## fit model using weighted sum of objectives
        function(w) {
          fitBy(MODEL,
            objective = function(...) sum(w * sapply(objective, function(obj) obj(...))),
            ...
          )
        }
      )
    },
    front <- apply(weights, 1, function(w) {
      fitBy(MODEL,
        objective = function(...) sum(w * sapply(objective, function(obj) obj(...))),
        ...
      )
    })
  ) ## switch parallel
  names(front) <- apply(weights, 1, paste, collapse = "_")
  front <- as.runlist(front)
  return(front)
}

#' @import utils
utils::globalVariables(c("cl"))
