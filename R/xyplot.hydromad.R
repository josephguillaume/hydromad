## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


plot.hydromad <-
  function(x, y, ...) {
    stop(
      "There is no 'plot' method for 'hydromad' objects.",
      "Try 'xyplot', or 'plot(fitted(...))'"
    )
  }

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

tsdiag.hydromad <- function(object, gof.lag, ...) {
  stats:::tsdiag.Arima(object$uh, gof.lag = gof.lag, ...)
}
