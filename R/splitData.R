## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Joseph Guillaume <josephguillaume@gmail.com>
##



#' Split data of a model instance
#'
#' Creates a runlist of model instances with separate data periods
#'
#'
#' @param object an object of class \code{hydromad}.
#' @param start.dates dates at which each model instance should start, passed
#' as \code{start} to \code{\link{window.zoo}}
#' @param end.dates dates at which each model instance should end, passed as
#' \code{end} to \code{\link{window.zoo}}
#' @param periods named list of start and end dates. Overrides
#' \code{start.dates} and \code{end.dates} if it is not \code{NULL}
#' @return A runlist of hydromad objects with data periods starting at
#' \code{start.dates} and ending at \code{end.dates}
#' @author Joseph Guillaume
#' @seealso \code{\link{crossValidate}}
#' @examples
#'
#' data(Cotter)
#' modx <- hydromad(Cotter,
#'   sma = "cwi", routing = "expuh",
#'   tau_s = c(2, 100), v_s = c(0, 1)
#' )
#' start.dates <- as.Date(c("1966-05-01", "1976-05-1"))
#' end.dates <- as.Date(c("1976-04-30", "1986-04-30"))
#' rl <- splitData(modx, start.dates, end.dates)
#' str(rl)
#' ## Fit the first period
#' fitx <- fitByOptim(rl[[1]])
#' ## Evaluate the resulting model on all periods
#' rl.first <- update(rl, newpars = coef(fitx))
#' @export
splitData <- function(object, start.dates = NULL, end.dates = NULL, periods = NULL) {
  ## Split data according to start.dates and end.dates
  if (!is.null(periods)) {
    start.dates <- sapply(periods, head, 1)
    end.dates <- sapply(periods, tail, 1)
  }
  cv.data <- mapply(window, start = start.dates, end = end.dates, MoreArgs = list(x = observed(object, all = TRUE, select = TRUE)), SIMPLIFY = FALSE)
  rl <- as.runlist(mapply(update, newdata = cv.data, MoreArgs = list(object = object), SIMPLIFY = FALSE))
  if (!is.null(names(periods))) {
    names(rl) <- names(periods)
  } else {
    names(rl) <- paste(start.dates, end.dates, sep = "_")
  }
  return(rl)
}
