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



#' Rolling cross-correlation at given lags.
#'
#' Rolling cross-correlation at given lags.  Can be useful to show how the
#' relationship between two time series changes over time, including out-by-one
#' timing errors.
#'
#' This is a fairly straightforward application of \code{\link{rollapply}} with
#' the \code{\link{ccf}} function. It may be better to do a time-varying
#' regression between the two series.
#'
#' @importFrom stats na.contiguous ccf
#' @importFrom lattice simpleTheme strip.default simpleKey
#' @importFrom latticeExtra xyplot.list
#' @importFrom zoo as.zoo rollapply
#'
#' @aliases rollccf xyplot.rollccf
#' @param DATA a named list, data frame, time series or zoo object containing
#' the two data series.
#' @param width a list or number specifying the width of window(s), in time
#' steps, in which to calculate cross correlations.
#' @param by temporal resolution: cross correlation is calculated in windows
#' spaced every \code{by} time steps.
#' @param lags,base.lag \code{lags} for which to calculate the cross
#' correlation. By default these are based on the overall maximum cross
#' correlation, \code{base.lag}.
#' @param rises if \code{TRUE}, compute the cross correlation with \emph{rises
#' in} streamflow. In this case the streamflow series must be named \code{"Q"}.
#' @param na.action function to handle missing data in each window (not the
#' whole series). This is only applied when the number of missing values is
#' less than \code{na.max.fraction}.
#'
#' Could be \code{na.exclude}.
#' @param na.max.fraction if the proportion of missing values in the moving
#' window exceeds this value, the corresponding result will be \code{NA}.
#' @param x,data \code{x} is an object produced by the \code{rollccf} function.
#' \code{data} is ignored.
#' @param with.data if \code{TRUE}, include the original data series in the
#' plot.
#' @param type,type.data drawing styles for the cross correlation series and
#' input data series. See \code{\link{panel.xyplot}}.
#' @param par.settings,layout,strip,ylim,xlab,as.table,...  passed to
#' \code{\link{xyplot}}.
#' @return \code{rollccf} returns a list of class \code{"rollccf"}, with
#' components: \item{rolls}{a list of time series of cross correlations.  One
#' element for each value of \code{width}. } \item{data}{time series input
#' data. } \item{lags, width, call}{ values of arguments used. }
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{ccf}}, \code{\link{rollapply}}
#' @keywords ts
#' @examples
#'
#' data(Canning)
#' foo <- rollccf(Canning)
#' xyplot(foo)
#' @export
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


#' @rdname rollccf
#' @export
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


#' @rdname rollccf
#' @export
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
