## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


event.xyplot <- function(x, ...) {
  UseMethod("event.xyplot")
}

event.xyplot.formula <-
  function(x, data = list(), events, ...,
           response = NULL, eFUN = sum,
           as.table = TRUE, layout = NULL,
           auto.key = NCOL(response) > 1,
           xlab = "event-aggregated covariate values",
           ylab = "event-aggregated response") {
    x <- as.formula(x)
    ## for zoo objects this splits up the columns into items:
    data <- as.list(data)
    ## define a special function to aggregate events, for use in formula
    e <- function(X, FUN = sum, ...) eventapply(X, events, FUN = FUN, ...)
    ## make it known to the formula:
    environment(x) <- environment()
    ## another convenience function for formula:
    first <- function(x) head(x, 1)
    ## extract covariates from formula
    tt <- terms(x, data = data, keep.order = TRUE)
    mf <- model.frame(tt, data, na.action = na.pass) ## na.omit stuffs up index
    if (is.null(response)) {
      response <- model.response(mf)
      if (is.null(response)) {
        stop("no response in formula, and 'response' argument missing")
      }
      mf <- mf[, -1, drop = FALSE]
    } else {
      response <- e(response, eFUN)
    }
    ystack <- stack(as.data.frame(response))
    foo <-
      xyplot.list(mf,
        data = ystack,
        FUN = function(VAR, ...) {
          xyplot(values ~ rep(coredata(VAR), length = NCOL(response) * NROW(response)),
            groups = ind, ...
          )
        },
        default.scales = list(x = list(relation = "free")),
        ...,
        as.table = as.table, layout = layout,
        auto.key = auto.key,
        xlab = xlab, ylab = ylab,
        y.same = NA
      ) ## allow specified 'ylim'
    foo$call <- sys.call(sys.parent())
    foo
  }
