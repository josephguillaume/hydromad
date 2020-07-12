## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


event.xyplot.hydromad.runlist <-
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
           panel.groups = function(x, y, ...) {
             panel.xyplot(x, y, ...)
             if (!requireNamespace("mgcv")) stop("package mgcv is required for event.xyplot")
             panel.smoother(x, y, y ~ s(x), method = "gam", ...)
           },
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
