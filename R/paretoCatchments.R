#' @importFrom utils write.csv
#' @importFrom stats aggregate
#' @importFrom reshape rename
#'
#' @export
paretoTimeAnalysis.data.frame <- function(stat, show.models = NA, objectives = "r.squared",
                                          pars, ...) {
  cat("\nCross-validation Pareto analysis\nWhich models cannot be rejected, due to dataset uncertainty/non-stationarity?\n")
  cat("\n== Eliminating Pareto-dominated models ==\n (worse than another model in all periods)\n")

  stat <- as.data.frame(stat)
  if (length(unique(stat$sim.period)) == 1) stop("Only one sim.period found, need at least two for this analysis to be make sense")
  if (is.null(stat$Catchment)) stat$Catchment <- "Somewhere"
  if (is.null(stat$Model.str)) stat$Model.str <- "A model"

  stat.split <- split(stat, stat$Catchment)
  stat.split <- lapply(stat.split, paretoTimeAnalysis_areModelsDominated, objectives = objectives)


  nondom.long <- do.call("rbind", lapply(stat.split, getDominatedLong))
  stopifnot(!any(is.na(nondom.long$dominated)))
  nondom.long <- nondom.long[!nondom.long$dominated, ]

  cat("\nHow many models are dominated (and therefore eliminated)\n")
  print(sapply(stat.split, function(res) table(res$dominated)), row.names = FALSE)

  cat("\nWhich model structures are non-dominated in each catchment?\nProportion of model instances that are non-dominated\n")
  ## model.ok <- do.call(rbind,lapply(stat.split,function(xx) sort(unique(stat$Model.str)) %in% xx$Model.str[!xx$dominated]))
  model.ok <- do.call(rbind, lapply(stat.split, function(xx) {
    xx$Model.str <- factor(xx$Model.str)
    signif(table(xx$Model.str[!xx$dominated]) / table(xx$Model.str), 3)
  }))
  print(model.ok, row.names = FALSE)

  if (!is.na(show.models)) {
    isd <- do.call(rbind, stat.split)
    ## FIXME isd$calib.period <- period.labels[as.character(isd$calib.period)]
    isd$dominated <- ifelse(isd$dominated, "Yes", "No")
    isd <- isd[
      with(isd, order(Catchment, dominated, Model.str, calib.period)),
      c(
        c("Catchment", "dominated", "Model.str", "calib.period"),
        setdiff(names(isd), c("Catchment", "dominated", "Model.str", "calib.period"))
      )
    ]
    if (isTRUE(show.models)) {
      cat("\nSummary of model realisations and whether they are dominated for each catchment\n")
      print(isd, row.names = FALSE)
    } else if (!is.na(show.models)) {
      cat("\nWriting csv of model realisations and whether they are dominated for each catchment\n")
      write.csv(isd, sprintf("%s_isdominated_models_catchments.csv", show.models), row.names = FALSE)
    }
  } else {
    cat("\nSpecify show.models=TRUE to show non-dominated and dominated models\nSpecify show.models=prefix to obtain csv of whether models are dominated\n")
  }

  cat("
== Performance across all periods ==
What is the range of non-dominated performance (RNDP) across all periods?
Is it large - is changing datasets causing problems?
")
  min.perf.all <- merge(
    aggregate(value ~ Catchment, data = nondom.long, min),
    aggregate(value ~ Catchment, data = nondom.long, max),
    by = c("Catchment"), suffixes = c(".min", ".max"), sort = FALSE
  )
  min.perf.all$RNDP <- min.perf.all$value.max - min.perf.all$value.min
  print(min.perf.all, row.names = FALSE)


  cat("
== Performance in each period ==
What is the RNDP in each period?
Is it low even though total RNDP is high? Why?
 Is there reason to believe the objective function is not comparable over time?
")
  min.perf.periods <- merge(
    aggregate(value ~ sim.period + Catchment, data = nondom.long, min),
    aggregate(value ~ sim.period + Catchment, data = nondom.long, max),
    by = c("Catchment", "sim.period"), suffixes = c(".min", ".max")
  )
  min.perf.periods <- min.perf.periods[with(
    min.perf.periods,
    order(Catchment, sim.period)
  ), ]
  min.perf.periods$RNDP <- min.perf.periods$value.max - min.perf.periods$value.min
  print(min.perf.periods, row.names = FALSE)

  cat("
== Worst non-dominated models in each period ==
Do any non-dominated models have unacceptable performance?
Which non-dominated model has the worst performance in each period? Why?
 Is it consistently the same dataset? Is there reason for that dataset to be problematic?
 Is it consistently the same model structure? Should another model structure have been selected?
 Is it consistently the same calibration objective function? Is it overfitting part of the record?
")
  for (cc in names(stat.split)) {
    cat("\n", cc, "\n")
    wmodel <- apply(stat.split[[cc]][!stat.split[[cc]]$dominated, as.character(unique(stat$sim.period))], 2, which.min)
    worst <- data.frame(
      sim.period = unique(stat$sim.period),
      worst.performance = apply(stat.split[[cc]][!stat.split[[cc]]$dominated, as.character(unique(stat$sim.period))], 2, min),
      stat.split[[cc]][!stat.split[[cc]]$dominated, ][wmodel, ]
    )
    worst$Catchment <- NULL # redundant column - used as title
    worst$dominated <- NULL # redundant column - subsetted to be false
    worst$objective <- NULL # redundant column - specified as argument
    print(worst, row.names = FALSE)
  }

  cat("
== Variation in inferred internal behaviour - Range of non-dominated parameters ==
Is the range of parameters of non-dominated models large for any single model structure?
Does the difference in performance correspond to different internal model behaviour?

")
  if (missing(pars)) {
    cat("argument 'pars' missing, skipping section\n")
  } else {
    nondom.idvars <- unique(subset(nondom.long, select = -c(sim.period, value)))
    coefs <- merge(nondom.idvars, pars)

    par.ranges <- merge(
      aggregate(value ~ variable + Model.str + Catchment, data = coefs, mean),
      aggregate(value ~ variable + Model.str + Catchment, data = coefs, min),
      by = c("Catchment", "Model.str", "variable"), suffixes = c(".mean", ".min"), sort = FALSE
    )
    par.ranges <- merge(
      par.ranges,
      aggregate(value ~ variable + Model.str + Catchment, data = coefs, max),
      by = c("Catchment", "Model.str", "variable"), sort = FALSE
    )
    par.ranges <- rename(par.ranges, c(value = "value.max"))
    par.ranges$range <- round(par.ranges$value.max - par.ranges$value.min, 3)
    par.ranges$"range.as.%.of.mean" <- round(par.ranges$range / par.ranges$value.mean * 100, 2)
    print(par.ranges, row.names = FALSE)
  }

  invisible(NULL)
}

getDominatedLong <- function(res.dom) {
  resm2 <- as.data.frame(res.dom)
  resm2 <- melt(resm2,
    id.vars = intersect(names(res.dom), c("Model.str", "Catchment", "calib.period", "Cal.objfn", "objective", "dominated")),
    variable_name = "sim.period"
  )
  ## FIXME: naming scheme of model realisations
  ## resm2$model <- paste(resm2$Model.str,resm2$calib.period,resm2$Cal.objfn)
  resm2
}

#' @export
paretoTimeAnalysis.matrix <- paretoTimeAnalysis.data.frame
