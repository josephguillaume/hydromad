## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## aspects:
## * aggregation groups (regular / events) & aggregation function
## * reference model
## * data transformation

buildTsObjective <-
  function(Q, groups = NULL, FUN = sum, ...,
           ref = NULL, boxcox = FALSE, start = NULL) {
    attributesQ <- attributes(Q)
    doaggr <- identity
    if (!is.null(groups)) {
      argsForFUN <- list(...)
      fullFUN <- function(...) {
        do.call(FUN, modifyList(list(...), argsForFUN))
      }
      doaggr <- function(x) {
        eventapply(x, groups, FUN = fullFUN)
      }
    }
    aggrQ <- doaggr(Q)
    aggrRef <- NULL
    if (!is.null(ref)) {
      aggrRef <- doaggr(ref)
    }
    if (!identical(boxcox, FALSE) && (length(aggrQ) > 1)) {
      coreaggrQ <- coredata(na.omit(aggrQ))
      if (is.null(start)) {
        start <-
          quantile(coreaggrQ[coreaggrQ > 0], 0.1, names = FALSE)
      }
      if (isTRUE(boxcox)) {
        lambda <- coef(powerTransform(coreaggrQ + start))
      } else {
        stopifnot(is.numeric(boxcox))
        lambda <- boxcox
      }
      function(Q, X, ...) {
        nseStat(aggrQ, doaggr(X),
          ref = aggrRef, ...,
          trans = function(x) bcPower(x + start, lambda)
        )
      }
    } else {
      function(Q, X, ...) {
        if (!missing(Q)) {
          if (!identical(attributes(Q), attributesQ)) {
            warning("'Q' has different attributes to that passed to buildTsObjectiveFun()")
          }
        }
        ## TODO: check that this 'Q' has same shape / index as original 'Q'?
        nseStat(aggrQ, doaggr(X), ref = aggrRef, ...)
      }
    }
  }

buildTsLikelihood <-
  function(Q, groups = NULL, FUN = sum, ...,
           boxcox = FALSE, start = NULL,
           distribution = dnorm, outliers = 0) {
    doaggr <- identity
    trans <- identity
    if (!is.null(groups)) {
      argsForFUN <- list(...)
      fullFUN <- function(...) {
        do.call(FUN, modifyList(list(...), argsForFUN))
      }
      doaggr <- function(x) {
        eventapply(x, groups, FUN = fullFUN)
      }
    }
    aggrQ <- doaggr(Q)
    if (!identical(boxcox, FALSE) && (length(aggrQ) > 1)) {
      coreaggrQ <- coredata(na.omit(aggrQ))
      if (is.null(start)) {
        start <-
          quantile(coreaggrQ[coreaggrQ > 0], 0.1, names = FALSE)
      }
      if (isTRUE(boxcox)) {
        lambda <- coef(powerTransform(coreaggrQ + start))
      } else {
        stopifnot(is.numeric(boxcox))
        lambda <- boxcox
      }
      trans <- function(x) {
        x[] <- bcPower(x + start, lambda)
        x
      }
    }
    function(Q, X) {
      resids <- trans(aggrQ) - trans(doaggr(X))
      logp <- distribution(resids, log = TRUE)
      if (outliers > 0) {
        logp <- tail(sort(logp), -outliers)
      }
      sum(logp)
    }
  }
