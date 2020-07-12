## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

fitStat <- function(...) {
  .Deprecated("nseStat")
  nseStat(...)
}

nseStat <-
  function(obs, mod, ref = NULL, ..., p = 2,
           trans = NULL, negatives.ok = FALSE,
           na.action = na.pass) {
    if (!is.vector(obs) || !is.vector(mod) ||
      (length(ref) > 1) && (!is.vector(ref))) {
      ## not plain vectors so assume time series and merge
      if (length(ref) > 1) {
        dat <- cbind(obs = obs, mod = mod, ref = ref)
      } else {
        dat <- cbind(obs = obs, mod = mod)
      }
      if (NROW(dat) <= 1) {
        warning("merged time series have no data; incompatible times?")
        return(NA_real_)
      }
      dat <- na.action(dat)
      if (NROW(dat) <= 1) {
        warning("time series have no data after 'na.action'")
        return(NA_real_)
      }
      obs <- dat[, "obs"]
      mod <- dat[, "mod"]
      if (length(ref) > 1) {
        ref <- dat[, "ref"]
      }
    } else {
      if (!identical(attributes(obs), attributes(mod))) {
        warning("attributes of 'obs' and 'mod' are not identical")
      }
    }
    obs <- coredata(obs)
    mod <- coredata(mod)
    stopifnot(length(obs) == length(mod))
    ## only use pairwise common data
    ok <- complete.cases(obs, mod)
    if (length(ref) > 1) {
      ref <- coredata(ref)
      stopifnot(length(ref) == length(mod))
      ok <- ok & !is.na(ref)
      ref <- ref[ok]
    }
    obs <- obs[ok]
    mod <- mod[ok]
    if (length(obs) == 0) {
      return(NA_real_)
    }
    ## negative values can cause errors with log()
    if (negatives.ok == FALSE) {
      obs <- pmax(obs, 0)
      mod <- pmax(mod, 0)
      if (!is.null(ref)) {
        ref <- pmax(ref, 0)
      }
    }
    ## transformation function
    if (!is.null(trans)) {
      if (is.character(trans)) {
        trans <- get(trans, mode = "function")
      }
      if (is.null(trans)) {
        trans <- identity
      }
      obs <- trans(obs)
      mod <- trans(mod)
      if (!is.null(ref)) {
        ref <- trans(ref)
      }
      ## check again, just in case 'trans' changed the length
      if (length(obs) != length(mod)) {
        stop("length(obs) != length(mod) after transformation")
      }
      ## check again for missing values introduced by 'trans'
      ok2 <- complete.cases(obs, mod)
      obs <- obs[ok2]
      mod <- mod[ok2]
      if (length(ref) > 1) {
        ref <- ref[ok2]
      }
    }
    ## default ref is the mean of transformed 'obs'
    if (is.null(ref)) {
      ref <- mean(obs)
    }
    ## calculate absolute error for model and reference
    ## and apply power p
    err <- abs(obs - mod)^p
    referr <- abs(obs - ref)^p
    1 - sum(err) / sum(referr)
  }
