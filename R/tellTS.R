# Call tell method for each column of ts.matrix, fun should return
#  long format data.frame: ts.id,variable,variables,...



#' Time series sensitivity analysis
#'
#' Generic method to evaluate sensitivity analysis methods in decoupled mode
#' (see \link[sensitivity]{tell}) on a time series of results, as produced for
#' example by \link{evalParsRollapply}. See examples for more details.
#'
# @importFrom parallel clusterEvalQ
#' @importFrom stats sd
#'
#' @aliases tellTS tellTS.default tellTS.sobol2002 tellTS.sobol2007
#' tellTS.morris
#' @param x A sensitivity analysis object for which results still need to be
#' provided. A typed list storing the state of the sensitivity study
#' (parameters, data, estimates), as returned by sensitivity analyses objects
#' constructors, such as morris, sobol2002, etc.
#' @param ts.matrix A matrix of model responses, with each row corresponding to
#' the parameter sets identified in \code{x} and each column a timestep for
#' which to evaluate sensitivity indices. If the matrix does not fit in memory,
#' \code{ff} matrices are also supported.
#' @param fun A function \code{fun(i,x)} which returns a \code{data.frame},
#' where \code{i} is the index of the time series and \code{x} is a complete
#' sensitivity analysis object, i.e. the result of calling \code{tell(x,y)}.
#' For sobol2002 and sobol2007, the default is to return the Total Sensitivity
#' Index (TSI) and its minimum and maximum bootstrap confidence interval.
#' @param indices (Optional) Which time indices to use. By default, sensitivity
#' indices are calculated for all columns of \code{ts.matrix}
#' @param parallel If "clusterApply", evaluate parameters in parallel using the
#' parallel or snow package. See \code{\link{hydromad_parallelisation}}
#' @param \dots Ignored by default implementation
#' @return a data.frame with columns produced by \code{fun}
#' @author Joseph Guillaume
#' @seealso \code{hydromad_sensitivity} for sensitivity analysis on
#' static rather than time series model output, \code{\link{evalParsRollapply}}
#' to obtain parameter sets
#' @references Herman, J. D., P. M. Reed, and T. Wagener. 2013. "Time-Varying
#' Sensitivity Analysis Clarifies the Effects of Watershed Model Formulation on
#' Model Behavior." Water Resources Research 49 (3): 1400-1414.
#' doi: \href{http://dx.doi.org/10.1002/wrcr.20124}{10.1002/wrcr.20124}
#' @keywords models
#' @examples
#'
#' library(sensitivity)
#'
#' ## Load data
#' data(Cotter)
#' obs <- Cotter[1:1000]
#'
#' ## Define rainfall-runoff model structure
#' model.str <- hydromad(obs,
#'   sma = "cwi", routing = "expuh",
#'   tau_q = c(0, 2), tau_s = c(2, 100), v_s = c(0, 1)
#' )
#'
#' ## Set the random seed to obtain replicable results
#' set.seed(19)
#'
#' ## Setup Morris Sensitivity analysis
#' incomplete <- morris(
#'   ## Names of factors/parameters
#'   factors = names(getFreeParsRanges(model.str)),
#'   ## Number of repetitions
#'   r = 4,
#'   ## List specifying design type and its parameters
#'   design = list(type = "oat", levels = 10, grid.jump = 2),
#'   ## Minimum value of each non-fixed parameter
#'   binf = sapply(getFreeParsRanges(model.str), min),
#'   ## Maximum value of each non-fixed parameter
#'   bsup = sapply(getFreeParsRanges(model.str), max)
#' )
#'
#'
#' # Calculate rolling time series of objective for each parameter set,
#' #   keeping results in memory
#' runs <- evalParsRollapply(incomplete$X, model.str,
#'   objective = ~ hmadstat("r.squared")(Q, X) /
#'     (2 - hmadstat("r.squared")(Q, X)),
#'   parallel = "none",
#'   filehash.name = NULL
#' )
#'
#' # Calculate Morris elementary effects for each timestep
#' sens <- tellTS(incomplete, runs, parallel = "none")
#'
#' ## Plot time series of mean absolute elementary effects
#' library(ggplot2)
#' qplot(x = ts.id, y = mu.star, colour = variable, data = sens, geom = "line")
#'
#' ## And of standard deviation of elementary effects
#' qplot(x = ts.id, y = sigma, colour = variable, data = sens, geom = "line")
#'
#' ######################################################################
#' \dontrun{
#' ## Sobol indices
#' ### Setup parallelisation
#' library(parallel)
#' hydromad.options("parallel" = list("evalParsRollapply" = "clusterApply"))
#' hydromad.options("parallel" = list("tellTS" = "clusterApply"))
#' cl <- makeCluster(3)
#' clusterEvalQ(cl, library(hydromad))
#'
#' ### Create sobol design
#' n <- 1000 ## Set number of samples desired
#' X1 <- parameterSets(getFreeParsRanges(model.str), n)
#' X2 <- parameterSets(getFreeParsRanges(model.str), n)
#' incomplete <- sobol2002(model = NULL, X1 = X1, X2 = X2, nboot = 100)
#'
#' ### Run model for each parameter set, storing result to disk
#' system.time(runs <- evalParsRollapply(incomplete$X, model.str,
#'   objective = hmadstat("r.squared"), parallel = "clusterApply"
#' ))
#'
#' ## Calculate sobol indices for each timestep with stored runs
#' system.time(TSI <- tellTS(incomplete, runs))
#'
#' ## Plot time series of sensitivities
#' library(ggplot2)
#' qplot(x = ts.id, y = original, colour = variable, data = TSI, geom = "line")
#'
#' ## Time series plot facetted by variable showing bootstrapped confidence interval
#' ggplot(data = TSI) +
#'   geom_ribbon(aes(
#'     x = ts.id, ymin = get("min. c.i."),
#'     ymax = get("max. c.i.")
#'   ), fill = "grey") +
#'   geom_line(aes(x = ts.id, y = original)) +
#'   facet_wrap(~variable)
#' }
#'
#' @export
tellTS <- function(x, ts.matrix, fun,
                   indices,
                   parallel = hydromad.getOption("parallel")[["tellTS"]],
                   ...) {
  UseMethod("tellTS")
}

#' @export
tellTS.default <- function(x, ts.matrix, fun,
                           indices,
                           parallel = hydromad.getOption("parallel")[["tellTS"]],
                           ...) {
  force(fun)
  if (is.null(parallel)) parallel <- "none"
  if (missing(indices)) indices <- 1:ncol(ts.matrix)
  cat(sprintf(
    "Calculating sensitivity on %d data points %s parallelisation\n", length(indices),
    ifelse(parallel != "none", "*without*", "with")
  ))
  switch(parallel,
    "clusterApply" = {
      parallel::clusterEvalQ(cl, requireNamespace("sensitivity"))
      parallel::clusterEvalQ(cl, requireNamespace("ff"))
      parallel::clusterExport(cl, c("ts.matrix", "x", "fun"), envir = environment())
      results <- parallel::parLapply(cl, indices, function(i) {
        y <- ts.matrix[, i]
        sensitivity::tell(x, y, ...)
        return(fun(i, x))
      })
    },
    {
      results <- lapply(indices, function(i) {
        y <- ts.matrix[, i]
        sensitivity::tell(x, y, ...)
        return(fun(i, x))
      })
    }
  )
  return(do.call(rbind, results))
}


#' @export
tellTS.sobol2002 <- function(x, ts.matrix, fun,
                             indices,
                             parallel = hydromad.getOption("parallel")[["tellTS"]],
                             ...) {
  if (missing(fun)) {
    fun <- function(i, x) {
      cbind(
        ts.id = i,
        variable = rownames(x$T),
        x$T[, c("original", "min. c.i.", "max. c.i.")]
      )
    }
  }

  tellTS.default(x, ts.matrix, fun,
    indices,
    parallel = hydromad.getOption("parallel")[["tellTS"]],
    ...
  )
}

#' @export
tellTS.sobol2007 <- tellTS.sobol2002

#' @export
tellTS.morris <- function(x, ts.matrix, fun,
                          indices,
                          parallel = hydromad.getOption("parallel")[["tellTS"]],
                          ...) {
  if (missing(fun)) {
    fun <- function(i, x) {
      res <- data.frame(
        ts.id = i,
        variable = colnames(x$ee),
        mu = apply(x$ee, 2, mean),
        mu.star = apply(x$ee, 2, function(x) mean(abs(x))),
        sigma = apply(x$ee, 2, sd)
      )
      rownames(res) <- NULL
      res
    }
  }

  tellTS.default(x, ts.matrix, fun,
    indices,
    parallel = hydromad.getOption("parallel")[["tellTS"]],
    ...
  )
}

#' @import utils
utils::globalVariables(c("cl"))
