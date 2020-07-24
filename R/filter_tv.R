## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


## time-varying recursive filter


#' Recursive filter with time-varying coefficient
#'
#' This is simply an AR(1) recursive filter with a time-varying recession rate.
#'
#'
#' @param x numeric vector or time series.
#' @param a numeric vector the same length as \code{x}, giving the filter
#' coefficient at each time step.
#' @param init value for \code{x[0]}.
#' @return a numeric vector or time series, like \code{x}.
#' @note If there are internal missing values, these are skipped over in the
#' calculation, maintaining the "state": the value of \code{y[i-1]} after any
#' missing values is the value from just before them. This behaviour is
#' different from \code{\link{filter}}, which drops the state back to 0.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{filter}}
#' @keywords ts
#' @examples
#'
#' ## The non-compiled function is this simple, if there are no NAs:
#' ftv2 <- function(x, a, init = 0) {
#'   y <- x
#'   y[1] <- x[1] + a[1] * init
#'   for (i in 2:length(x)) {
#'     y[i] <- x[i] + a[i] * y[i - 1]
#'   }
#'   return(y)
#' }
#' ## make a sine wave filter
#' a <- sin(pi * seq(0, 3 * pi, length = 100)) * 0.2 + 0.9
#' plot.ts(a, ylim = c(0, 1.2))
#' ## response to a unit impluse
#' x <- c(1, rep(0, 99))
#' y <- filter_tv(x, a)
#' plot.ts(y, log = "y")
#' stopifnot(isTRUE(all.equal(y, ftv2(x, a))))
#' ## treatment of missing values
#' x[15:20] <- NA
#' plot.ts(filter_tv(x, a), log = "y")
#' @useDynLib hydromad ar1_tv
#' @export
filter_tv <-
  function(x, a, init = 0) {
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(a))
    stopifnot(is.numeric(init))
    stopifnot(length(x) == length(a))
    ## skip over missing values (maintaining the state y[i-1])
    bad <- !is.finite(x) | !is.finite(a)
    x[bad] <- 0
    a[bad] <- 1
    if (hydromad.getOption("pure.R.code")) {
      y <- x
      y[1] <- x[1] + a[1] * init
      for (i in 2:length(x)) {
        y[i] <- x[i] + a[i] * y[i - 1]
      }
    } else {
      y <- x
      y[] <- .C(ar1_tv,
        as.double(x),
        as.double(a),
        as.integer(length(x)),
        as.double(init),
        out = double(length(x)),
        PACKAGE = "hydromad"
      )$out
    }
    ## re-insert missing values
    y[bad] <- NA
    y
  }

## recursive filter skimming over NAs (treated as zeros)
## TODO: or could maintain state using na.exclude(), naresid()
filter_ok <-
  function(x, filter, method = "recursive", ...) {
    bad <- !is.finite(x)
    x[bad] <- 0
    y <- filter(x, filter = filter, method = method, ...)
    ## re-insert missing values
    y[bad] <- NA
    y
  }
