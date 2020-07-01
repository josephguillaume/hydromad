## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


## time-varying recursive filter
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
        DUP = FALSE, PACKAGE = "hydromad"
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
