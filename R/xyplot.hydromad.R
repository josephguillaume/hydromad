## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

#' Plot methods for Hydromad model objects
#'
#' Plot methods...
#'
#' @importFrom lattice qqmath xyplot make.groups
#' @importFrom stats tsdiag ppoints
#' @importFrom latticeExtra as.layer xyplot.list
#' @importFrom zoo coredata
#'
#' @name xyplot.hydromad
#' @aliases xyplot.hydromad.runlist qqmath.hydromad plot.hydromad
#' tsdiag.hydromad
#' @param x an object of class \code{hydromad}.
#' @param data ignored.
#' @param scales Placeholder
#' @param type Placeholder
#' @param layout Placeholder
#' @param auto.key Placeholder
#' @param \dots further arguments passed on to \code{\link{xyplot.zoo}} or
#' \code{\link{qqmath}}.
#' @param feasible.bounds if \code{TRUE}, then ensemble simulation bounds are
#' extracted and plotted. This only works if a \emph{feasible set} has been
#' specified using \code{\link{defineFeasibleSet}} or the \code{update} method.
#' Note that the meaning depends on what value of \code{glue.quantiles} was
#' specified to those methods: it might be the overall simulation bounds, or
#' some GLUE-like quantile values.
#' @param col.bounds,border,alpha.bounds graphical parameters of the ensemble
#' simulation bounds if \code{feasible.bounds = TRUE}.
#' @param all passed to \code{fitted()} and \code{observed()}.
#' @param with.P to include the input rainfall series in the plot.
#' @param type.P plot type for rainfall, passed to \code{\link{panel.xyplot}}.
#' @param superpose to overlay observed and modelled time series in one panel.
#' @param f.value,tails.n arguments to \code{\link{panel.qqmath}}.
#' @param object,gof.lag passed to the \code{arima} method of
#' \code{\link{tsdiag}}.
#' @param y Placeholder for plot.hydromad
#' @return the trellis functions return a trellis object.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{hydromad.object}, \code{\link{xyplot}},
#' \code{\link{xyplot.ts}}, \code{\link{xyplot.list}}
#' @keywords hplot ts
#' @examples
#'
#' data(Canning)
#' cannCal <- window(Canning, start = "1978-01-01", end = "1982-12-31")
#' mod <-
#'   hydromad(cannCal,
#'     sma = "cwi", tw = 162, f = 2, l = 300,
#'     t_ref = 0, scale = 0.000284,
#'     routing = "expuh", tau_s = 4.3, delay = 1, warmup = 200
#'   )
#'
#' xyplot(mod, with.P = TRUE)
#'
#' c(
#'   streamflow = xyplot(mod),
#'   residuals = xyplot(residuals(mod, type = "h")),
#'   layout = c(1, 2), y.same = TRUE
#' )
#'
#' xyplot(residuals(mod)) +
#'   latticeExtra::layer(panel.tskernel(..., width = 90, c = 2, col = 1)) +
#'   latticeExtra::layer(panel.tskernel(..., width = 180, c = 2, col = 1, lwd = 2)) +
#'   latticeExtra::layer(panel.tskernel(..., width = 360, c = 2, lwd = 2))
#'
#' qqmath(mod,
#'   scales = list(y = list(log = TRUE)), distribution = qnorm,
#'   type = c("g", "l")
#' )
#' @export
plot.hydromad <-
  function(x, y, ...) {
    stop(
      "There is no 'plot' method for 'hydromad' objects.",
      "Try 'xyplot', or 'plot(fitted(...))'"
    )
  }


#' @rdname xyplot.hydromad
#' @export
xyplot.hydromad <-
  function(x, data = NULL, ..., scales = list(),
           feasible.bounds = FALSE,
           col.bounds = "grey80", border = "grey60", alpha.bounds = 1,
           all = FALSE, superpose = TRUE,
           with.P = FALSE, type = "l",
           type.P = c("h", if ("g" %in% type) "g"),
           layout = c(1, NA)) {
    stopifnot(is.null(data))

    if (isValidModel(x)) {
      tsdat <- cbind(
        observed = observed(x, all = all),
        modelled = fitted(x, all = all)
      )
    } else {
      tsdat <- observed(x, all = all)
    }
    foo <- xyplot(tsdat, ...,
      scales = scales, superpose = superpose,
      type = type
    )
    if (feasible.bounds) {
      bounds <- fitted(x, all = all, feasible.bounds = TRUE)
      if (NCOL(bounds) > 2) { ## show lowest and highest quantiles
        bounds <- bounds[, c(1, NCOL(bounds))]
      }
      ## make a whole plot, passing scales, rather than just a layer
      ## because the y scale may be log in which case the data must be transformed.
      foo <- foo +
        as.layer(xyplot(bounds, ...,
          scales = scales, superpose = TRUE, type = type,
          col = col.bounds, alpha = alpha.bounds, border = border,
          panel = function(x, y, ...) {
            x2 <- matrix(x, ncol = 2)
            y2 <- matrix(y, ncol = 2)
            panel.ribbon(zoo(y2, x2[, 1]), ...)
          }
        ),
        under = TRUE
        )
    }

    if (with.P) {
      ## never want rainfall to be on a log scale:
      scales$y$log <- FALSE
      rainPlot <-
        xyplot(observed(x, select = "P", all = all), ...,
          scales = scales, superpose = superpose,
          type = type.P
        )
      foo <- c(
        streamflow = foo, rainfall = rainPlot,
        x.same = NA, y.same = NA, layout = layout
      )
    }
    foo$call <- sys.call(sys.parent())
    foo
  }


#' @rdname xyplot.hydromad
#' @export
xyplot.hydromad.runlist <-
  function(x, data = NULL, ..., scales = list(),
           all = FALSE, superpose = FALSE,
           with.P = FALSE, type = "l",
           type.P = c("h", if ("g" %in% type) "g"),
           layout = c(1, NA)) {
    stopifnot(is.null(data))
    if (superpose) {
      ## fitted models superposed, but rainfall still juxtaposed.
      ## include observed series from item 1 (assuming all are the same!)
      tsdat <- cbind(
        observed = observed(x[[1]], all = all),
        fitted(x, all = all)
      )
      foo <- xyplot(tsdat, ...,
        superpose = superpose, scales = scales,
        type = type, layout = layout
      )
    } else {
      ## fitted models juxtaposed, each with observed flow superposed
      foo <- xyplot.list(x, ...,
        all = all, scales = scales,
        superpose = TRUE, with.P = FALSE,
        type = type, layout = layout
      )
    }
    if (with.P) {
      ## never want rainfall to be on a log scale:
      scales$y$log <- FALSE
      rainPlot <-
        xyplot(observed(x[[1]], select = "P", all = all), ...,
          scales = scales, superpose = superpose,
          type = type.P
        )
      foo <- c(foo,
        rainfall = rainPlot,
        x.same = TRUE, y.same = NA, layout = layout
      )
      if (superpose) {
        rownames(foo)[1] <- "streamflow"
      }
    }
    foo$call <- sys.call(sys.parent())
    foo
  }


#' @rdname xyplot.hydromad
#' @export
qqmath.hydromad <-
  function(x, data = NULL, ..., all = FALSE, type = "l",
           auto.key = list(lines = TRUE, points = FALSE),
           f.value = ppoints(100), tails.n = 100) {
    stopifnot(is.null(data))
    tsdat <- cbind(
      obs = observed(x, all = all),
      mod = fitted(x, all = all)
    )
    ## keep only common (corresponding) values
    coredata(tsdat)[complete.cases(tsdat) == FALSE, ] <- NA
    dat <- make.groups(observed = tsdat[, "obs"], modelled = tsdat[, "mod"])
    foo <- qqmath(~data,
      groups = which, data = dat,
      f.value = f.value, tails.n = tails.n,
      auto.key = auto.key, type = type, ...
    )
    foo$call <- sys.call(sys.parent())
    foo
  }


#' @rdname xyplot.hydromad
#' @export
tsdiag.hydromad <- function(object, gof.lag, ...) {
  tsdiag.Arima <- getS3method("tsdiag", "Arima")
  tsdiag.Arima(object$uh, gof.lag = gof.lag, ...)
}
