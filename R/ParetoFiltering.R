
paretoFilter <- function(x, ...) {
  d <- ncol(x)
  n <- nrow(x)
  is.optimal <- rep(TRUE, n)
  for (i in 1:(n - 1)) {
    for (j in i:n) {
      if (i != j && (is.optimal[i] || is.optimal[j])) {
        xi <- x[i, ]
        xj <- x[j, ]
        if (all(xi <= xj) && any(xi < xj)) {
          is.optimal[j] <- FALSE
        }
        else if (all(xj <= xi) && any(xj < xi)) {
          is.optimal[i] <- FALSE
        }
      }
    }
  }
  is.optimal
}

plotPCNSE <- function(res, objectives = "r.squared", return.data = FALSE) {
  if (!require("ggplot2")) stop("package ggplot2 is required for plotPCNSE")
  res <- as.data.frame(res)
  stopifnot(!"Catchment" %in% names(res) | length(unique(res$Catchment)) == 1)
  resm <- melt(res,
    measure.vars = objectives,
    id.vars = intersect(names(res), c("Model.str", "Catchment", "calib.period", "sim.period", "Cal.objfn")),
    variable_name = "objective"
  )
  res2 <- cast(resm, ... ~ sim.period)
  sim.periods <- unique(resm$sim.period)

  res2$dominated <- !paretoFilter(-as.matrix(res2[, sim.periods]))
  resm2 <- as.data.frame(res2)
  resm2 <- melt(resm2, id.vars = intersect(names(res2), c("Model.str", "Catchment", "calib.period", "Cal.objfn", "objective", "dominated")), variable_name = "sim.period")

  ## Niceties
  ## resm2$sim.period <- ordered(resm2$sim.period,
  ##                             levels=intersect(c("70","80","90","00"),unique(resm2$sim.period))
  ##                             )
  ## FIXME resm2$sim.period <- ordered(period.labels[as.character(resm2$sim.period)])
  resm2$sim.period <- ordered(resm2$sim.period)
  ## resm2$model <- paste(resm2$Model.str,resm2$Catchment,resm2$calib.period)
  ## FIXME resm2$model <- paste(resm2$Model.str,resm2$calib.period,resm2$Cal.objfn)
  resm2$model <- paste(resm2$Model.str, resm2$calib.period)
  resm2$dominated <- ifelse(resm2$dominated, "Yes", "No")

  ## browser()
  if (return.data) {
    return(resm2)
  }

  ggplot(resm2) +
    geom_line(aes(x = sim.period, y = value, group = model, col = dominated, linetype = Model.str)) +
    ## geom_line(aes(x=sim.period,y=value,group=model),col="grey")+
    geom_point(aes(x = sim.period, y = value, group = model, col = dominated)) +
    geom_text(aes(x = sim.period, y = value, label = model),
      size = 3, hjust = -0.2,
      data = subset(resm2, sim.period == max(resm2$sim.period))
    ) +
    geom_text(aes(x = sim.period, y = value, label = model),
      size = 3, hjust = 1.2,
      data = subset(resm2, sim.period == min(resm2$sim.period))
    ) +
    scale_y_continuous(name = "Performance - NSE", breaks = seq(0.4, 1, by = 0.1)) +
    scale_x_discrete(name = "Simulation period") +
    scale_colour_discrete(name = "Dominated?") +
    scale_linetype_discrete(name = "Model structure") +
    coord_cartesian(ylim = c(0.4, 1)) ## FIXME
}
