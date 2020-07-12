## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


# rollccf2 <-
#    function(DATA = list(Q =, P=), ...)
# {
## TODO: better: use dlm

## regression with basis splines over time
## -- approach suggested by Trevor Hastie
## and posted by Tim Hesterberg on R-help 2007-01-23
#    library("splines")
#    timebasis <- bs(time, df = df)
#    fit <- lm(y ~ timebasis * x, DATA)
#    slcoef <- coef(fit)[-(1:(df+1))]
#    slope <- slcoef[1] + timebasis %*% slcoef[-1]
## TODO
#    slope
# }

rollccf <-
  function(DATA = data.frame(Q = , P = ),
           width = list(365, 90),
           by = 28,
           lags = base.lag + c(0, 1, -1),
           base.lag = estimateDelay(DATA, rises = rises, plot = FALSE),
           rises = FALSE,
           na.action = na.contiguous,
           na.max.fraction = 1 / 3) {
    DATA <- as.zoo(DATA)
    stopifnot(NCOL(DATA) >= 2)
    if (rises) {
      if ("Q" %in% colnames(DATA)) {
        DATA[, "Q"] <- pmax(diff(DATA[, "Q"], na.pad = TRUE), 0)
      } else {
        stop("Give an item 'Q' to use flow rises.")
      }
    }
    if ("Q" %in% colnames(DATA)) {
      ## make 'Q' the first column, so that positive lags are physical
      whichQ <- which(colnames(DATA) == "Q")
      DATA <- DATA[, c(whichQ, (1:NCOL(DATA))[-whichQ])]
    }
    DATA <- DATA[, 1:2]
    width <- as.list(width)
    if (is.null(names(width))) {
      names(width) <- paste("width", unlist(width))
    }
    obj <- list()
    obj$rolls <-
      lapply(
        width,
        function(wid) {
          tmp <- rollapply(DATA,
            width = wid, by = by,
            by.column = FALSE,
            FUN = ccfForLags, lags = lags,
            na.action = na.action,
            na.max.fraction = na.max.fraction
          )
          tmp
        }
      )
    obj$data <- DATA
    obj$lags <- lags
    obj$width <- width
    obj$call <- match.call()
    class(obj) <- c("rollccf", class(obj))
    obj
  }

xyplot.rollccf <-
  function(x, data = NULL, ...,
           with.data = TRUE,
           type = list(c("h", "b")),
           type.data = "l",
           par.settings = simpleTheme(pch = ".", cex = 2),
           layout = c(1, length(x$rolls) + with.data * 2),
           strip = strip.default,
           ylim = c(0, 1), xlab = NULL, as.table = TRUE) {
    rollplot <- xyplot.list(x$rolls,
      x.same = TRUE, y.same = NA,
      type = type, ylim = ylim,
      superpose = TRUE,
      par.settings = par.settings, xlab = xlab,
      as.table = as.table,
      key = simpleKey(paste("lag", x$lags),
        lines = TRUE, points = FALSE, columns = 3
      )
    )
    if (with.data) {
      rollplot <- c(rollplot,
        xyplot(x$data,
          type = type.data,
          par.settings = par.settings
        ),
        x.same = TRUE
      )
    }
    rollplot <- update(rollplot, strip = strip, ...)
    rownames(rollplot) <- c(names(x$rolls), if (with.data) colnames(x$data))
    rollplot$layout <- layout
    rollplot$call <- sys.call(sys.parent())
    return(rollplot)
  }

ccfForLags <- function(DATA, lags = 0,
                       na.action = na.contiguous,
                       na.max.fraction = 1 / 3) {
  vals <- rep(NA, length(lags))
  names(vals) <- paste("lag", lags)
  if ((sum(complete.cases(DATA)) <= 10) ||
    (mean(complete.cases(DATA)) < na.max.fraction)) {
    return(vals)
  }
  ans <- ccf(DATA[, 1], DATA[, 2],
    lag.max = max(abs(lags)),
    plot = FALSE, na.action = na.action
  )
  vals[] <- drop((ans[lags])$acf)
  vals
}
