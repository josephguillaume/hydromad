#' Analysis using Pareto-filtering of model performance across simulation
#' periods
#'
#' Identify dominated model realisations, that are inferior to another model in
#' all simulation periods, optionally in multiple catchments.  Summary results
#' are then produced to help assess how performance of model realisations and
#' model structures varies across simulation periods.
#'
#' @name paretoTimeAnalysis
#' @aliases paretoTimeAnalysis_performance
#' paretoTimeAnalysis_areModelsDominated paretoTimeAnalysis.data.frame
#' paretoTimeAnalysis.matrix
#' @param \dots arguments to \code{paretoTimeAnalysis.crossvalidation} or
#' \code{paretoTimeAnalysis.data.frame}
#' @param rl a crossvalidation \code{runlist}: a list of fitted model objects
#' produced by \code{\link{crossValidate}}
#' @param show.models If \code{TRUE}, print out a table of models and an
#' indication of whether they are dominated, as produced by
#' \code{paretoTimeAnalysis_areModelsDominated}. If not \code{NA}, write out
#' the table to "show.models_isdominated_models_catchments.csv")
#' @param objectives Vector of column names containing performance measures
#' used to determine whether models are dominated across time periods. We
#' assume higher values are better. Values should be transformed prior to use.
#' @param qoi Quantities of interest to calculate, as interpreted by
#' \code{\link{objFunVal}}. Defaults are given as examples: 90th percentile
#' runoff (a scalar prediction) and R Squared using log transformed data (a
#' performance statistic).
# @param stat data.frame of results, including the column \code{sim.period}
# and the columns named in \code{objectives}. The following columns may be
# included as id variables:
# \code{Model.str},\code{Catchment},\code{calib.period},\code{Cal.objfn}
# @param pars long format data.frame of parameter values for each model with
# columns \code{variable} and \code{value}.  The following columns may be
# included as id variables:
# \code{Model.str},\code{Catchment},\code{calib.period},\code{Cal.objfn}.  If
# \code{pars} is missing, results that require it are skipped.
#' @return For \code{paretoTimeAnalysis}, no return value. Used for its side
#' effect of producing text. Optionally writes csv files (see the argument
#' \code{show.models}).
#'
#' \code{paretoTimeAnalysis_areModelsDominated} produces a wide-format
#' data.frame with id variable columns, a column for each \code{sim.period}
#' value, and a column \code{dominated} indicating whether another model is
#' better in all simulation periods.
#' @author Joseph Guillaume
#' @seealso \code{\link{crossValidate}} for an example use of
#' \code{paretoTimeAnalysis.crossvalidation}
#' @references Placeholder
#' @keywords models
#' @examples
#'
#' ## Dataset consisting of results for two simulation periods,
#' ##  obtained by calibration in the same periods with different
#' ##  model structures.
#' data(YeAl97)
#'
#' ## For one catchment, produce a table indicating whether models defined by
#' ## their calib.period and Model.str are dominated according to the objective E
#' paretoTimeAnalysis_areModelsDominated(subset(YeAl97, Catchment == "Salmon"), objectives = "E")
#'
#' ## For all catchments, performance analysis
#' paretoTimeAnalysis(YeAl97, objectives = "E")
#' @export
paretoTimeAnalysis <- function(...) UseMethod("paretoTimeAnalysis")

#' @rdname paretoTimeAnalysis
#' @export
paretoTimeAnalysis.crossvalidation <- function(rl,
                                               show.models = NA,
                                               objectives = "r.squared",
                                               qoi = list(
                                                 q90.in.units.of.runoff = function(Q, X, ...) quantile(X, 0.9, na.rm = TRUE),
                                                 r.sq.log = hmadstat("r.sq.log")
                                               ), ...) {
  stopifnot(inherits(rl, "crossvalidation"))
  stat <- summary(rl, ...)
  pars <- coef(rl)
  pars <- melt(unique(subset(pars, select = -sim.period)),
    na.rm = TRUE,
    id.vars = intersect(names(pars), c("calib.period", "Model.str", "Cal.objfn", "Catchment"))
  )
  paretoTimeAnalysis.data.frame(stat, show.models = show.models, objectives = objectives, pars = pars)

  cat("
== Performance with other statistics ==
Use the set of non-dominated models as an ensemble to predict a quantity of interest
Is the model unacceptable in any period? Is the uncertainty too large?
")

  ## TODO: avoid recomputation of dominated models cf paretoCatchments
  stat <- as.data.frame(stat)
  if (is.null(stat$Catchment)) stat$Catchment <- "Somewhere"
  stat.split <- split(stat, stat$Catchment)
  stat.split <- do.call(rbind, lapply(stat.split, paretoTimeAnalysis_areModelsDominated, objectives = objectives))
  id.vars <- intersect(names(stat), c(
    "Model.str", "Catchment",
    "calib.period", "Cal.objfn"
  ))
  is.nondominated <- apply(stat[, id.vars], 1, paste, collapse = "_") %in% apply(stat.split[!stat.split$dominated, id.vars], 1, paste, collapse = "_")
  stopifnot(length(is.nondominated) == length(rl))

  ## Calculate qois for each model
  p <- lapply(rl[is.nondominated], objFunVal, objective = qoi)
  ## Convert from list of lists to data.frame
  p <- as.data.frame(do.call(rbind, lapply(p, unlist)))
  ## Extract simulation and calibration period names
  p$sim.period <- sub("_cal.*", "", rownames(p))
  p$calib.period <- sub(".*_cal", "", rownames(p))
  ## Melt qois
  pm <- melt(p, id.vars = intersect(names(p), c("Model.str", "Catchment", "sim.period", "calib.period", "Cal.objfn", "objective", "dominated")))
  ## Cast and aggregate showing min, max and range of each qoi
  p2 <- cast(pm, variable + sim.period ~ .,
    fun.aggregate = function(x) c(min = min(x), max = max(x), range = diff(range(x)))
  )
  print(p2, row.names = FALSE)

  invisible(NULL)
} ## paretoTimeAnalysis
