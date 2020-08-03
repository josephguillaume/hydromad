#' Pareto filter
#'
#' Filter a matrix of values to identify pareto-optimal solutions
#'
#'
#' @param x matrix of values with each row representing a solution
#' @param \dots ignored
#' @return Those values in 'x' which are not dominated by any other solution.
#' @author From the \code{mco} package.  Heike Trautmann <email:
#' trautmann@@statistik.uni-dortmund.de>, Detlef Steuer <email:
#' steuer@@hsu-hamburg.de> and Olaf Mersmann <email:
#' olafm@@statistik.uni-dortmund.de>
#' @seealso \code{\link{paretoTimeAnalysis_areModelsDominated}} which uses this
#' function to evaluate model performance across time-periods
#' @keywords models
#' @examples
#'
#' ## Performance measures from 4 models in the Salmon catchment,
#' ##  see YeAl97
#' mat <- matrix(c(
#'   0.865, 0.892, -0.847, 0.795,
#'   0.774, 0.905, 0.819, 0.930
#' ), nrow = 4)
#' mat
#' ## Identify dominated rows of the matrix, interpreting
#' ##   higher values to be better
#' ## TRUE: The 2nd and 4th rows are pareto-optimal/non-dominated
#' ## FALSE: The 1st and 3rd rows are both inferior to the 2nd row,
#' ##   and are therefore dominated
#' paretoFilter(-mat)
#' @export paretoFilter
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



#' Annotated parallel coordinates plot of model performance across periods
#'
#' Plot performance of model realisations, identifying non-dominated models
#'
#' @importFrom reshape melt cast
#'
# @importFrom ggplot2 geom_point geom_text scale_y_continuous scale_x_discrete
# scale_linetype_discrete scale_colour_discrete coord_cartesian ggplot
# aes geom_line
#'
#' @param res data.frame of results, including the column \code{sim.period} and
#' the columns named in \code{objectives}. At least one of the following column
#' names should be included as id variables:
#' \code{Model.str},\code{Catchment},\code{calib.period},,\code{Cal.objfn}
#' Other columns will be ignored.
#' @param objectives Vector of column names containing performance measures. We
#' assume higher values are better. Values should be transformed prior to use.
#' @param return.data return data.frame used in call to \code{ggplot}.
#' Facilitates custom plotting.
#' @return \code{ggplot} object, which can be plotted. Or data.frame if
#' \code{return.data} is \code{TRUE}.
#' @author Joseph Guillaume
#' @seealso \code{areModelsDominated} for raw data,
#' \code{paretoCatchments} for further analysis
#' @keywords models
#' @examples
#'
#' data(YeAl97)
#' plotPCNSE(subset(YeAl97, Catchment == "Salmon"), objectives = "E")
#' @export
plotPCNSE <- function(res, objectives = "r.squared", return.data = FALSE) {
  if (!requireNamespace("ggplot2")) stop("package ggplot2 is required for plotPCNSE")
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

  ggplot2::ggplot(resm2) +
    ggplot2::geom_line(ggplot2::aes(x = sim.period, y = value, group = model, col = dominated, linetype = Model.str)) +
    ## geom_line(aes(x=sim.period,y=value,group=model),col="grey")+
    ggplot2::geom_point(ggplot2::aes(x = sim.period, y = value, group = model, col = dominated)) +
    ggplot2::geom_text(ggplot2::aes(x = sim.period, y = value, label = model),
      size = 3, hjust = -0.2,
      data = subset(resm2, sim.period == max(resm2$sim.period))
    ) +
    ggplot2::geom_text(ggplot2::aes(x = sim.period, y = value, label = model),
      size = 3, hjust = 1.2,
      data = subset(resm2, sim.period == min(resm2$sim.period))
    ) +
    ggplot2::scale_y_continuous(name = "Performance - NSE", breaks = seq(0.4, 1, by = 0.1)) +
    ggplot2::scale_x_discrete(name = "Simulation period") +
    ggplot2::scale_colour_discrete(name = "Dominated?") +
    ggplot2::scale_linetype_discrete(name = "Model structure") +
    ggplot2::coord_cartesian(ylim = c(0.4, 1)) ## FIXME
}


#' @import utils
utils::globalVariables(c("sim.period", "value", "model", "dominated", "Model.str"))
