## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


## This plots fitted vs observed only.
## To plot residuals, just call xyplot(residuals())


#' Plot results from a set of model runs
#'
#' Plot results from a set of model runs using Lattice graphics.
#'
#' @importFrom lattice qqmath xyplot make.groups
#' @importFrom latticeExtra xyplot.list
#'
#'
#' @aliases xyplot.runlist qqmath.runlist
#' @param x a \code{runlist} object, which is a list of fitted model objects.
#' @param data ignored.
#' @param \dots further arguments are passed on to plotting functions.
#' @param all passed to \code{fitted()} and \code{observed()}.
#' @param superpose to overlay all model result series in one panel.
#' @param x.same,y.same passed to \code{\link{xyplot.list}}.
#' series.
#' @param residuals to plot the residual series rather than fitted and observed
#' series.
#' @param f.value,tails.n arguments to \code{\link{panel.qqmath}}.
#' @param layout Placeholder
#' @param auto.key Placeholder
#' @param type Placeholder
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{runlist}}, \code{\link{xyplot.ts}},
#' \code{\link{xyplot}}, \code{\link{qqmath}}
#' @keywords utilities
#' @examples
#'
#' data(HydroTestData)
#' mod1 <- hydromad(HydroTestData,
#'   sma = "scalar",
#'   routing = "expuh", tau_s = 10
#' )
#' mod2 <- update(mod1, tau_s = 20, tau_q = 5, v_s = 0.5)
#' mod3 <- update(mod2, loss = 0.5)
#' mods <-
#'   runlist(
#'     `single store` = mod1,
#'     `two stores` = mod2,
#'     loss = mod3
#'   )
#' xyplot(mods, superpose = TRUE)
#' xyplot(mods, scales = list(y = list(log = TRUE)))
#' @export
xyplot.runlist <-
  function(x, data = NULL, ...,
           all = FALSE, superpose = FALSE,
           x.same = TRUE, y.same = NA, layout = c(1, NA)) {
    if (superpose) {
      ## include observed series from item 1 (assuming all are the same!)
      tsdat <- cbind(
        observed = observed(x[[1]], all = all),
        fitted(x, all = all)
      )
      foo <- xyplot(tsdat, ..., superpose = TRUE, layout = layout)
    } else {
      foo <- xyplot.list(x, ...,
        x.same = x.same, y.same = y.same, layout = layout
      )
    }
    foo$call <- sys.call(sys.parent())
    foo
  }

## Handles either fitted vs observed, or residuals.


#' @rdname xyplot.runlist
#' @export
qqmath.runlist <-
  function(x, data = NULL, ..., all = FALSE,
           residuals = FALSE, superpose = FALSE,
           f.value = ppoints(100), tails.n = 100, type = "l",
           auto.key = list(lines = TRUE, points = FALSE)) {
    if (superpose) {
      if (residuals) {
        tsdat <- residuals(x, all = all)
      } else {
        ## include observed series from item 1 (assuming all are the same!)
        tsdat <- cbind(
          observed = observed(x[[1]], all = all),
          fitted(x, all = all)
        )
      }
      dat <- do.call("make.groups", as.data.frame(tsdat))
      foo <- qqmath(~data,
        groups = which, data = dat, ...,
        f.value = f.value, tails.n = tails.n,
        type = type, auto.key = auto.key
      )
    } else {
      if (residuals) {
        FUN <- function(x, ...) qqmath(~ residuals(x), ...)
      } else {
        FUN <- qqmath
      }
      foo <- xyplot.list(x, ...,
        FUN = FUN,
        f.value = f.value, tails.n = tails.n,
        type = type, auto.key = auto.key
      )
    }
    foo$call <- sys.call(sys.parent())
    foo
  }
