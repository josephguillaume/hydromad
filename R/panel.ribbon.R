##
## Copyright (c) 2008 Felix Andrews <felix@nfrac.org>
##

panel.ribbon <- function(...) {
  UseMethod("panel.ribbon")
}

## Plot the area between 2 series as a filled polygon.
## With groups, acts like panel.superpose, but with polygon style settings.
panel.ribbon.default <-
  function(x, y, y2, groups = NULL,
           col = if (is.null(groups)) plot.polygon$col else superpose.polygon$col,
           border = if (is.null(groups)) plot.polygon$border else superpose.polygon$border,
           lty = if (is.null(groups)) plot.polygon$lty else superpose.polygon$lty,
           lwd = if (is.null(groups)) plot.polygon$lwd else superpose.polygon$lwd,
           alpha = if (is.null(groups)) plot.polygon$alpha else superpose.polygon$alpha,
           ..., col.line = border, fill, panel.groups = panel.ribbon) {
    plot.polygon <- trellis.par.get("plot.polygon")
    superpose.polygon <- trellis.par.get("superpose.polygon")
    x <- as.numeric(x)
    y1 <- as.numeric(y)
    y2 <- as.numeric(y2)
    ylim <- current.panel.limits()$y
    infi <- is.infinite(y1)
    y1[infi] <- ifelse(y1[infi] > 0, max(ylim), min(ylim))
    infi <- is.infinite(y2)
    y2[infi] <- ifelse(y2[infi] > 0, max(ylim), min(ylim))
    if (length(x) == 0) {
      return()
    }
    if (!is.null(groups)) {
      ## NOTE superpose does not handle 'border' argument, so pass it as col.line
      ## NOTE superpose only subsets 'x' and 'y' so need to rename arguments:
      panel.superpose(y1, y2, ...,
        groups = groups, trueX = x,
        panel.groups = function(x, y, ..., trueX) {
          panel.ribbon.default(x = trueX, y = x, y2 = y, ...)
        },
        col = col, col.line = border, lty = lty, lwd = lwd,
        alpha = alpha
      )
    } else {
      if (all(is.na(col)) && !missing(col.line)) {
        col <- col.line
      }
      if (is.unsorted(x, na.rm = TRUE)) {
        stop("'x' should be ordered (increasing)")
      }
      ## need to split up the series into chunks without any missing values
      ## (because NAs split the polygon)
      xyy <- data.frame(x = x, y1 = y1, y2 = y2)
      ok <- complete.cases(xyy)
      runs <- rle(ok)
      ## assign unique values to each chunk, and NAs between (dropped by 'split')
      runs$values[runs$values == TRUE] <- seq_len(sum(runs$values))
      runs$values[runs$values == FALSE] <- NA
      ## expand into long format
      chunks <- inverse.rle(runs)
      lapply(split(xyy, chunks), function(xyyi, ...) {
        x <- xyyi$x
        y1 <- xyyi$y1
        y2 <- xyyi$y2
        ## join up the polygon
        xx <- c(x, rev(x))
        yy <- c(y1, rev(y2))
        ## we need to catch the 'fill' argument from panel.superpose, otherwise over-rides 'col'
        panel.polygon(xx, yy, alpha = alpha, col = col, border = col.line, lty = lty, lwd = lwd, ...)
      }, ...)
    }
  }

panel.ribbon.ts <- function(y, y2 = NULL, ...) {
  ## allow one 'ts' argument with 2 columns
  if (is.null(y2)) {
    if (NCOL(y) > 1) {
      y2 <- y[, 2]
      y <- y[, 1]
    } else {
      stop("y2 missing and not enough columns in y")
    }
  }
  panel.ribbon(x = as.vector(time(y)), y, y2, ...)
}

panel.ribbon.zoo <- function(y, y2 = NULL, ...) {
  ## allow one 'zoo' argument with 2 columns
  if (is.null(y2)) {
    if (NCOL(y) > 1) {
      y2 <- y[, 2]
      y <- y[, 1]
    } else {
      stop("y2 missing and not enough columns in y")
    }
  }
  panel.ribbon(x = index(y), coredata(y), coredata(y2), ...)
}
