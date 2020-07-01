## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


## This plots fitted vs observed only.
## To plot residuals, just call xyplot(residuals())
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
