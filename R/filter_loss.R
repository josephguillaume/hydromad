## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


## recursive filter with a constant loss
filter_loss <-
  function(x, a, loss, init = 0) {
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(a))
    stopifnot(is.numeric(loss))
    stopifnot(is.numeric(init))
    stopifnot(length(x) > length(a) + 1)
    stopifnot(length(a) > 0)
    init <- rep(init, length = length(a))
    xAttrs <- attributes(x)
    ## skip over missing values (maintaining the state y[i-1])
    bad <- !is.finite(x)
    x[bad] <- 0
    ## define function for a single (univariate) time series
    filter_loss_onecol <- function(x) {
      if (hydromad.getOption("pure.R.code")) {
        nn <- length(a)
        y <- c(init, x * 0)
        x <- c(init * 0, x)
        for (i in (nn + 1):length(x)) {
          y[i] <- x[i] + a * y[i - (1:nn)]
          y[i] <- max(0, y[i] - loss)
        }
        y <- y[-(1:nn)]
      } else {
        nn <- length(a)
        yi <- c(init, x * 0)
        x <- c(init * 0, x)
        y <- .C(filter_constloss,
          as.double(x),
          as.integer(length(x)),
          as.double(a),
          as.integer(length(a)),
          as.double(loss),
          out = as.double(yi),
          DUP = FALSE, PACKAGE = "hydromad"
        )$out
        y <- y[-(1:nn)]
      }
      y
    }
    if (NCOL(x) == 1) {
      y <- filter_loss_onecol(x)
    }
    if (NCOL(x) > 1) {
      y <- x
      for (k in 1:NCOL(x)) {
        y[, k] <- filter_loss_onecol(x[, k])
      }
    }
    ## make it a time series object again
    attributes(y) <- xAttrs
    ## re-insert missing values
    y[bad] <- NA
    y
  }
