## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


observed <- function(object, ...) {
  UseMethod("observed")
}

observed.default <- function(object, ...) {
  fitted(object, ...) + residuals(object, ...)
}



autocorrTime <-
  function(x, order = c(2, 1), fallTo = 0.001) {
    stopifnot(length(order) == 2)
    stopifnot((0 < fallTo) && (fallTo < 1))
    mod <- suppressWarnings(arima(x,
      order = c(order[1], 0, order[2]),
      include.mean = FALSE
    ))
    pars <- coef(mod)
    ar <- pars[seq_len(order[1])]
    ma <- pars[seq_len(order[2]) + order[1]]
    ans <- match(TRUE, ARMAtoMA(ar, ma, lag.max = NROW(x) / 2) < fallTo)
    if (is.na(ans)) {
      ans <- NROW(x) / 2
    }
    ans <- max(5, ans)
    ans
  }

lagbind <-
  function(x, lags, all = TRUE, dframe = FALSE) {
    bind.fn <- if (all) ts.union else ts.intersect
    names(lags) <- paste("t", ifelse(lags < 0, "-", "+"),
      abs(lags),
      sep = ""
    )
    do.call(bind.fn, c(
      lapply(lags, function(k) lag(x, k)),
      list(dframe = dframe)
    ))
  }


## TODO: replace shiftWindow use with lagts when it appears in zoo
shiftWindow <-
  function(x, delay, and.lag = FALSE, fill = NA) {
    ## this is faster than a full ts.intersect - just shift by delay time steps
    ## a positive delay shifts the window back in time,
    ## a negative delay shifts the window forward in time
    if (delay == 0) {
      return(x)
    }
    attribs <- attributes(x)
    x <- coredata(x)
    if (NCOL(x) == 1) {
      res <- if (delay > 0) {
        c(rep(fill, delay), x[-(length(x) - seq(delay) + 1)])
      } else {
        c(x[-seq(abs(delay))], rep(fill, abs(delay)))
      }
    } else {
      res <- if (delay > 0) {
        x[c(rep(NA, delay), seq(1, NCOL(x) - delay)), ]
      } else {
        x[c(seq(abs(delay) + 1, NCOL(x)), rep(NA, abs(delay))), ]
      }
      res[1:delay, ] <- fill
    }
    attributes(res) <- attribs
    if (and.lag == FALSE) {
      ## values are already implicitly lagged above;
      ## the following keeps values in sync with the time index:
      res <- lag(res, delay)
    }
    res
  }

# stripWarmup <-
#    function(x, warmup)
# {
#    UseMethod("stripWarmup")
# }

stripWarmup <-
  function(x, warmup) {
    if (length(x) == 0) {
      return(x)
    }
    if (warmup == 0) {
      return(x)
    }
    stopifnot(warmup >= 0)
    stopifnot(warmup < NROW(x) + 1)
    newstart <- time(x)[warmup + 1]
    window(x, start = newstart)
  }

# stripWarmup.default <-
#    function(x, warmup)
# {
#    if (length(x) == 0) return(x)
#    if (warmup == 0) return(x)
#    stopifnot(warmup >= 0)
#    stopifnot(warmup < NROW(x)+1)
#    if (!is.ts(x)) x <- as.ts(x)
#    freq <- frequency(x)
#    t_warm <- tsp(x)[1] + (warmup) / freq
#    t_end <- tsp(x)[2]
#    ## following equivalent to return(window(x, start=t_warm)) but faster(?)
#    if (NCOL(x) == 1)
#        x <- x[-(1:warmup)]
#    else 	x <- x[-(1:warmup),]
#    ans <- ts(x, t_warm, t_end, freq)
# }

numericSummary <- function(x) {
  x <- as.data.frame(x)
  foo <- sapply(x, function(z) {
    ans <- c(
      fivenum(z, na.rm = TRUE),
      mean(z, na.rm = TRUE)
    )
    names(ans) <- c(
      "Min.", "1st Qu.", "Median",
      "3rd Qu.", "Max.", "Mean"
    )
    ## place mean after median
    ans <- ans[c(1:3, 6, 4:5)]
    nas <- sum(is.na(z))
    ans["NA's"] <- nas
    ans
  })
  # if (all(foo["NA's",] == 0)) ...
  foo
}
