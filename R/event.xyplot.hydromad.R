## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##



#' Visualise systematic model errors against covariates
#'
#' Visualise systematic model errors against covariates.
#'
#' @importFrom lattice panel.xyplot trellis.par.get current.panel.limits
#' panel.superpose panel.polygon
#' @importFrom latticeExtra panel.smoother
#'
#' @name event.xyplot.hydromad
#' @aliases event.xyplot.hydromad.runlist
#' @param x a \code{hydromad} or \code{hydromad.runlist} object.
#' @param events event sequence produced by \code{\link{eventseq}}, or a vector
#' defining continguous groups on the specified variables.
#' @param formula formula defining the covariates to plot, as passed to the
#' formula method of \code{\link{event.xyplot}}. It may refer to any of the
#' variables in the model data frame (\code{observed(x, select = TRUE)});
#' additionally it may refer to the value \code{U}, for the effective rainfall
#' series derived from the model.  Finally \code{\link{julian}} may be referred
#' to (from data time index).
#' @param extract a function to apply to \code{x} to extract the response
#' variable; by default this is \code{residuals} but could be e.g.
#' \code{fitted} or \code{function(x) residuals(x, boxcox = TRUE)}.
#' @param with.U to include modelled effective rainfall \code{U} as a
#' covariate.
#' @param \dots further arguments passed to \code{\link{event.xyplot}} and on
#' to \code{\link{xyplot}} and the panel function.
#' @param panel,panel.groups,abline,pch,ylab passed to \code{\link{xyplot}}.
#' @param data ignored.
#' @return this function returns a trellis object which can be \code{plot}ted.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{event.xyplot}}, \code{\link{xyplot}},
#' \code{\link{eventapply}}
#' @keywords hplot ts
#' @examples
#'
#' data(Cotter)
#' x <- Cotter[1:1000, ]
#' mod <- hydromad(x,
#'   sma = "scalar",
#'   routing = "armax", rfit = list("sriv", order = c(2, 1))
#' )
#' ev <- eventseq(x$P, thresh = 3, inthresh = 1, indur = 5)
#' event.xyplot(mod, events = ev)
#' event.xyplot(mod,
#'   events = ev,
#'   extract = function(x) residuals(x, boxcox = TRUE)
#' )
#'
#' foo <- event.xyplot(mod,
#'   events = ev,
#'   ~ sqrt(e(P, max)) + sqrt(e(rollmean(lag(P, -1), 20, align = "left"), first))
#' )
#' dimnames(foo)[[1]] <- c("sqrt. peak rain (mm/day)", "mean 20-day ante. rain")
#' foo
#' @export
event.xyplot.hydromad <-
  function(x, events,
           formula =
             ~ log2(e(Q, mean) + .01) +
               log2(e(lag(Q, -2), first) + .01) +
               log2(e(U, max) + .01) +
               e(E, mean),
           extract = residuals,
           with.U = TRUE,
           ...,
           panel = panel.superpose,
           panel.groups = panel.groups.funs,
           abline = list(h = 0), pch = ".",
           ylab = "residual flow sums in event windows (mm)",
           data = NULL) {
    if (!is.null(data)) warning("'data' ignored")
    if (inherits(x, "runlist")) {
      data <- observed(x[[1]], select = TRUE)
    } else {
      data <- observed(x, select = TRUE)
    }
    if (with.U) {
      data$U <- fitted(x, U = TRUE)
    }
    data$julian <- julian(time(data))
    response <- extract(x)
    foo <-
      event.xyplot(formula,
        data = data,
        events = events,
        response = response,
        ...,
        panel = panel, panel.groups = panel.groups,
        abline = abline, pch = pch,
        ylab = ylab
      )
    foo$call <- sys.call(sys.parent())
    foo
  }

#' @rdname event.xyplot.hydromad
#' @export
event.xyplot.hydromad.runlist <- event.xyplot.hydromad

panel.groups.funs <- function(x, y, ...) {
  panel.xyplot(x, y, ...)
  if (!requireNamespace("mgcv")) stop("package mgcv is required for event.xyplot")
  panel.smoother(x, y, y ~ s(x), method = "gam", ...)
}
