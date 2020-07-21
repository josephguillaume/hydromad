## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


# Deprecated
fitStat <- function(...) {
  .Deprecated("nseStat")
  nseStat(...)
}



#' Generalisation of Nash-Sutcliffe Efficiency (R Squared)
#'
#' Generalisation of Nash-Sutcliffe Efficiency (R Squared) as a fit statistic
#' for time series.
#'
#' The result is, after transformation of variables,
#'
#' \deqn{1 - sum(abs(obs-mod)^p) / sum(abs(obs-ref)^p)}
#'
#' A perfect fit gives a value of 1 and a fit equivalent to the reference model
#' gives a value of 0. Values less than 0 are worse than the reference model.
#'
#' If the arguments \code{obs}, \code{mod} or \code{ref} are not plain vectors,
#' \code{nseStat} will attempt to merge them together, so that corresponding
#' time steps are compared to each other even if the time windows are not
#' equal.
#'
#' @aliases nseStat fitStat
#' @param obs observed data vector.
#' @param mod model-predicted data vector corresponding to \code{obs}.
#' @param ref reference model predictions corresponding to \code{obs}. If
#' \code{NULL}, \code{ref} is taken as the mean of \code{obs} after applying
#' any transformation (\code{trans}).
#' @param ...  ignored.
#' @param p power to apply to absolute residuals (\code{abs(obs - mod)} and
#' \code{abs(obs - ref)}. The default, \code{p = 2} corresponds to the
#' Nash-Sutcliffe Efficiency (NSE). Setting \code{p = 1} does not square the
#' residuals and is sometimes called Normalised sum of Absolute Errors (NAE).
#' @param trans a function to apply to each data series before calculating the
#' fit statistic.
#' @param negatives.ok if \code{FALSE}, the default case, all values in
#' \code{obs}, \code{mod} and \code{ref} are constrained to be non-negative;
#' i.e. negative values are replaced with zero.
#' @param na.action a function to apply to the time series, which is expected
#' to fill in or remove missing values (note, this is optional).
#' @return a single numeric value.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad.stats}}, \code{\link{buildTsObjective}},
#' \code{\link{objFunVal}}, \code{\link{summary.hydromad}}
#' @keywords ts
#' @examples
#'
#' ## generate some data -- it is autocorrelated so the use of these
#' ## stats is somewhat problematic!
#' set.seed(0)
#' U <- ts(pmax(0, rgamma(200, shape = 0.1, scale = 20) - 5))
#' ## simulate error as multiplicative uniform random
#' Ue <- U * runif(200, min = 0.5, max = 1.5)
#' ## and resample 10 percent of time steps
#' ii <- sample(seq_along(U), 20)
#' Ue[ii] <- rev(U[ii])
#' ## apply recursive filter
#' Q <- filter(U, 0.7, method = "r")
#' X <- filter(Ue, 0.75, method = "r")
#'
#' ## convert to 'zoo' objects with Date index
#' Q <- zoo(Q, as.Date("2000-01-01") + 1:200)
#' X <- zoo(X, time(Q))
#'
#' xyplot(merge(Q, X), superpose = TRUE)
#'
#' nseStat(Q, X)
#'
#' nseStat(Q, X, trans = sqrt)
#'
#' nseStat(Q, X, trans = function(x) log(x + 1))
#'
#' ## use absolute residuals rather than squared residuals
#' nseStat(Q, X, p = 1)
#'
#' ## use a different reference model (one-step-ahead forecast)
#' nseStat(Q, X, ref = lag(Q, -1))
#'
#' ## reference as seasonal averages rather than overall average
#' nseStat(Q, X, ref = ave(Q, months(time(Q))))
#'
#' ## see how the reference model performs in terms of R Squared
#' nseStat(Q, ave(Q, months(time(Q))))
#' @export
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
