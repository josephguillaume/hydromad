## Expects columns:
## Catchment calib.period sim.period Model.str     + objectives + other
## Objectives need to be specified (default r.squared)
## Assume objectives higher value is better. Values should be transformed prior to use
paretoTimeAnalysis_areModelsDominated <- function(stat, objectives = "r.squared") {
  stat <- as.data.frame(stat)

  if (!all(objectives %in% names(stat))) {
    stop(sprintf("objectives not recognised: %s",
      paste(setdiff(objectives, names(stat))),
      collapse = ", "
    ))
  }

  id.vars <- intersect(names(stat), c("Model.str", "Catchment", "calib.period", "sim.period", "Cal.objfn"))
  stat <- stat[, c(id.vars, objectives)]

  statm <- melt(stat,
    id.vars = id.vars,
    measure.vars = objectives,
    variable_name = "objective"
  )
  stat2 <- cast(statm, ... ~ sim.period)
  sim.periods <- unique(statm$sim.period)

  stat2$dominated <- !paretoFilter(-as.matrix(stat2[, as.character(sim.periods)]))
  stat2
}
