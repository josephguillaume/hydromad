##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


#' Identify discrete events from time series and apply functions to them.
#'
#' Identify discrete events from time series and apply functions to them.
#'
#' @importFrom zoo is.regular rollmax
#'
#' @name events
#' @aliases eventseq eventapply eventinfo findThresh
#' @param x,X a \code{\link{ts}} or \code{\link{zoo}} object.  May be
#' multivariate, i.e. have multiple columns.
#' @param thresh threshold value: the data must be strictly above this level
#' (or, if \code{below = TRUE}, below this level) to define an event. If
#' \code{x} is a matrix (with multiple columns), \code{thresh} is allowed to be
#' a vector with one value per column. In this case, events continue while any
#' series exceeds the threshold. Finally, \code{thresh} can be a full matrix or
#' vector of the same dimension as \code{x}. Any missing values are treated as
#' below the threshold.
#' @param mingap the minimum number of time steps that can separate events. Any
#' inter-event durations shorter than this will be subsumed into the
#' surrounding event.
#' @param mindur the minimum number of time steps in each event window. Any
#' events whose duration is shorter than this will be skipped.
#' @param extend a number of time steps to include on the end of each event
#' window, before accounting for \code{inthresh}.  If this causes events to run
#' into following events they will be merged.
#' @param inthresh \code{inthresh} gives a second threshold value to define
#' when events stop: after an event is initiated by exceeding \code{thresh}, it
#' continues until the second (lower) threshold \code{inthresh} is reached (for
#' at least \code{indur} steps).  If missing or \code{NULL}, it defaults to
#' \code{thresh}.  As with \code{thresh}, \code{inthresh} is allowed to be a
#' vector with values corresponding to the columns of \code{inx}.
#' @param inx optionally, a different series may be given to determine event
#' termination with \code{inthresh}. E.g. an input series may define event
#' starts and a response series may define when they end.  If missing or
#' \code{NULL}, defaults to \code{x}.  As with \code{x}, this is allowed to be
#' a matrix.
#' @param indur the series must remain below \code{inthresh} for this many time
#' steps in order to terminate an event.
#' @param below reverses the definition of events, to identify events where the
#' data falls below \code{thresh} (and/or \code{inthresh}). Setting this to
#' \code{TRUE} is equivalent to negating \code{x} and \code{thresh}, etc.
#' @param all to include the periods between events as additional levels (i.e.
#' as events in themselves). These inter-events will be assigned negative
#' levels, while the actual events will have positive levels. Otherwise -- in
#' the default case -- these inter-event periods are left as \code{NA}.
#' @param continue set \code{continue = TRUE} to extend each event until the
#' following event; otherwise there are gaps between events (the default
#' behaviour). This has no effect if \code{all = TRUE}.
#' @param n number of events to be identified: if this is given (and > 0), then
#' a value of \code{thresh} is estimated such that about this many events are
#' returned. A warning is given if the actual number of events is not within 5
#' percent of \code{n}. Note that there may be multiple possible threshold
#' levels giving the same number of events!
#' @param events a \code{factor}-like object, perhaps \code{\link{ts}} or
#' \code{\link{zoo}}. Its levels (unique values) define events, and \code{NA}
#' values are ignored. If \code{events} is a factor or vector, its values are
#' assumed to correspond to \code{X}; otherwise, it is assumed to be a time
#' series object and is merged (\code{cbind}ed) with \code{X}.
#' @param FUN, a function to apply to the values in each event. In the default
#' case of \code{by.column = TRUE}, \code{FUN} will be passed a vector of
#' values from one column of \code{X}. When \code{by.column = FALSE}, which
#' only makes sense when \code{X} is a matrix, \code{FUN} will be passed a
#' matrix. The dots (\code{\dots{}}) are passed on too.
# @param list() a function to apply to the values in each event. In the
# default case of \code{by.column = TRUE}, \code{FUN} will be passed a vector
# of values from one column of \code{X}. When \code{by.column = FALSE}, which
# only makes sense when \code{X} is a matrix, \code{FUN} will be passed a
# matrix. The dots (\code{\dots{}}) are passed on too.
#' @param by.column a function to apply to the values in each event. In the
#' default case of \code{by.column = TRUE}, \code{FUN} will be passed a vector
#' of values from one column of \code{X}. When \code{by.column = FALSE}, which
#' only makes sense when \code{X} is a matrix, \code{FUN} will be passed a
#' matrix. The dots (\code{\dots{}}) are passed on too.
#' @param simplify if \code{FALSE}, the result will be returned as a list with
#' one (named) element for each event, rather than a time series like object.
#' This case allows \code{FUN} to return a complex object or vectors of
#' variable lengths.
#' @param TIMING defines how to construct the time index of the result. Should
#' the time corresponding to an event be taken from the \code{time()} of its
#' start, middle, or end?
#' @param within Placeholder
#' @param trace Placeholder
#' @param optimize.tol Placeholder
#' @param ... Placeholder
#' @return
#'
#' \code{eventseq} returns a zoo object, with core data consisting of an
#' ordered \code{\link{factor}}, representing the identified events, and the
#' same time index as \code{x}. Periods between events are left as \code{NA},
#' unless \code{all = TRUE} in which case they are treated as separate events.
#' The returned object stores \code{thresh} as an attribute.
#'
#' \code{eventapply} returns a \code{zoo} object (an irregular time series in
#' this case), with the value returned from \code{FUN} applied to each discrete
#' event in \code{X}.
#'
#' \code{eventinfo} returns a \code{data.frame} with columns \describe{
#' \item{list("Time")}{ time that the event started (from \code{time(X)}).  }
#' \item{list("Month, Year")}{ month and year (as integers) of the mid-point of
#' the event.  } \item{list("Value")}{ result of \code{FUN} applied to the
#' event.  } \item{list("Duration")}{ length of the event in time steps / data
#' points.  } \item{list("PreDuration")}{ number of time steps since the last
#' event ended.  } }
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{cut.Date}}, \code{\link{tapply}},
#' \code{\link{rollapply}}, \code{\link{aggregate.zoo}},
#' \code{\link{panel.xblocks}}, \code{clusters} in the \pkg{evd} package.
#' @keywords ts utilities
#' @examples
#'
#' data(Queanbeyan)
#' ## wet period
#' x <- window(Queanbeyan, start = "1974-01-01", end = "1976-12-01")
#'
#' evp <- eventseq(x$P, thresh = 5, inthresh = 1, indur = 4, continue = TRUE)
#' evq <- eventseq(x$Q, thresh = 2, indur = 4, mingap = 5)
#'
#' nlevels(evp) ## number of events
#' nlevels(evq)
#' str(evq)
#' table(coredata(evq))
#' eventapply(x$Q, evq, FUN = sum)
#' eventapply(x, evq, FUN = mean)
#' eventinfo(x$Q, evq)
#'
#' evplot <- xyplot(x) +
#'   latticeExtra::layer(panel.xblocks(evq, block.y = 0, vjust = 1, col = 1)) +
#'   latticeExtra::layer(panel.xblocks(evp, col = c("grey90", "grey80"), border = "grey80"))
#'
#' evplot
#'
#' update(evplot,
#'   type = "s",
#'   xlim = as.Date(c("1990-07-01", "1990-08-31"))
#' ) +
#'   latticeExtra::layer(panel.abline(h = c(5, 1), lty = 2), packets = 1)
#'
#' ## example of requesting a threshold giving about 'n' events
#' set.seed(0)
#' ee <- eventseq(rnorm(100), n = 10, mingap = 2)
#' nlevels(ee)
#' attr(ee, "thresh")
#'
#' ##
#' ## example of classifying events based on hydro properties
#' ##
#'
#' data(Queanbeyan)
#' x <- window(Queanbeyan, start = "1974-01-01", end = "1976-12-01")
#' e <- eventseq(x$P, thresh = 5, inthresh = 1, indur = 4, continue = TRUE)
#'
#' ## classify events based on max flow
#' qclass <- cut(
#'   ave(coredata(x$Q), coredata(e), FUN = max),
#'   c(0, 0.5, 1, Inf)
#' )
#' qclass <- zoo(qclass, time(x))
#'
#' ## Classify events based on antecedent flow
#' x <- merge(x, Q1 = lag(x$Q, -1), all = c(TRUE, FALSE))
#' head1 <- function(z) z[1]
#' q1class <- cut(
#'   ave(coredata(x$Q1), coredata(e), FUN = head1),
#'   c(0, 0.2, 0.3, Inf)
#' )
#' q1class <- zoo(q1class, time(x))
#'
#' ## combined classification
#' combin <- factor(paste("M", unclass(qclass), "_A", unclass(q1class), sep = ""))
#' combin <- zoo(combin, time(x))
#'
#' ## check results
#' head(data.frame(x, event = unclass(e), qclass, q1class, combin), 50)
#'
#' ## number of events in each class
#' each.e <- !duplicated(e)
#' table(coredata(combin[each.e]))
#' @export
findThresh <-
  function(x, thresh = NA, ## ignored
           n, within = (n %/% 20) + 1,
           mingap = 1, mindur = 1,
           below = FALSE, ...,
           trace = FALSE, optimize.tol = 0.1) {
    stopifnot(n > 0)
    n <- round(n)
    stopifnot(n * mindur + n * mingap <= NROW(x))
    x <- coredata(x)
    if (below) {
      x <- -x
    }
    ## return (difference from 'n' of) number of events for 'thresh'
    nDiffForThresh <- function(thresh) {
      ev <- eventseq(x,
        thresh = thresh, mingap = mingap,
        mindur = mindur, ...
      )
      newn <- nlevels(ev)
      if (trace) {
        message(sprintf("thresh = %.3f, n = %d", thresh, newn))
      }
      abs(newn - n)
    }

    ## do a quick run with a rough guess for the range
    rng <- quantile(x, (1 - (n * mindur / NROW(x)))^c(1, 3),
      na.rm = TRUE, names = FALSE
    )
    if (trace) {
      message("initial range guess: ", toString(signif(rng, 4)))
    }
    res <- optimize(nDiffForThresh,
      interval = rng,
      tol = optimize.tol
    )
    if (res$objective > within) {
      if (trace) {
        message("initial run off by ", res$objective)
      }
      ## expand possible range; only exclude minimum and maximum values
      rng2 <- quantile(x, range(ppoints(length(x))),
        na.rm = TRUE, names = FALSE
      )
      res2 <- optimize(nDiffForThresh,
        interval = rng2,
        tol = optimize.tol
      )
      if (res2$objective < res$objective) {
        res <- res2
      }
    }
    if (res$objective > within) {
      warning(
        "actual number of events differs from target (", n, ") by ",
        res$objective
      )
    }
    thresh <- res$minimum
    if (below) {
      thresh <- -thresh
    }
    return(thresh)
  }


#' @rdname events
#' @export
eventseq <-
  function(x, thresh = 0, mingap = 1, mindur = 1, extend = 0,
           inthresh = thresh, inx = x, indur = 1,
           below = FALSE, all = FALSE, continue = FALSE,
           n = NA) {
    if (!is.na(n) && (n > 0)) {
      if (thresh > 0) {
        warning("'thresh' argument ignored since 'n' given")
      }
      ccall <- match.call()
      ccall[[1]] <- quote(findThresh)
      thresh <- eval.parent(ccall)
    }
    ## check for simple case:
    inthreshGiven <-
      (((!missing(inthresh) || !missing(inx)) &&
        (!is.null(inthresh) && !is.null(inx))) ||
        (indur > 1))
    if (is.null(inthresh)) inthresh <- thresh
    if (is.null(inx)) inx <- x
    if (below) {
      x <- -x
      thresh <- -thresh
      inthesh <- -inthresh
      inx <- -inx
    }
    stopifnot(NROW(inx) == NROW(x))
    ## assume NAs are below threshold
    x[is.na(x)] <- -Inf
    inx[is.na(inx)] <- -Inf
    ## find runs above threshold
    if (is.matrix(x)) {
      ## (runs continue while any series/column is above thresh)
      ## threshold can be a vector, one for each column
      threshmat <- thresh
      if (length(thresh) > 1) {
        threshmat <- matrix(thresh, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
      }
      isover <- (rowSums(coredata(x) > threshmat) > 0)
    } else {
      isover <- (coredata(x) > thresh)
    }
    ## extend all events by 'extend' steps
    if (extend > 0) {
      ends <- c(FALSE, diff(isover) == -1)
      topad <- ends
      for (i in seq_len(extend - 1)) {
        topad <- topad | c(rep(FALSE, i), head(ends, -i))
      }
      isover[topad] <- TRUE
      #        uruns <- rle(isover)
      #        ii <- uruns$values == TRUE
      #        ## ignore any event at end
      #        ii[length(ii)] <- FALSE
      #        ## figure out how many extra steps can be taken from gap:
      #        ii.next <- which(ii)+1
      #        gap <- uruns$lengths[ii.next]
      #        okpad <- pmin(extend, gap)
      #        uruns$lengths[ii] <- uruns$lengths[ii] + okpad
      #        uruns$lengths[ii.next] <- uruns$lengths[ii.next] - okpad
      #        isover <- inverse.rle(uruns)
    }
    ## ensure that runs extend until 'inthresh' is met (inx <= inthresh)
    ## and remains below inthresh for at least 'indur' steps
    if (inthreshGiven) {
      if (indur > 1) { ## set inx to rolling maximum of last 'indur' steps
        inx <- rollmax(inx, indur, align = "right", fill = NA)
        inx <- na.locf(inx, fromLast = TRUE)
      }
      if (is.matrix(inx)) {
        ## threshold can be a vector, one for each column
        inthreshmat <- inthresh
        if (length(inthresh) > 1) {
          inthreshmat <- matrix(inthresh, nrow = nrow(inx), ncol = ncol(inx), byrow = TRUE)
        }
        stillover <- (rowSums(coredata(inx) > inthreshmat) > 0)
      } else {
        stillover <- (coredata(inx) > inthresh)
      }
      ends <- c(FALSE, diff(isover) == -1)
      ## from each of ends, while stillover, set isover = TRUE
      for (i in which(ends & stillover)) {
        j <- i
        while (j <= length(isover)) {
          if (!stillover[j]) break
          if (isover[j]) break ## ran into another event
          isover[j] <- TRUE
          j <- j + 1
        }
      }
    }
    ## run length encoding
    uruns <- rle(isover)
    ## find drops (between runs) whose length is less than 'mingap'
    if (mingap > 1) {
      nongap <- with(uruns, values == FALSE & lengths < mingap)
      ## ignore any gaps at start and end
      nongap[c(1, length(nongap))] <- FALSE
      if (any(nongap)) {
        ## set short gaps to be part of the surrounding cluster
        uruns$values[nongap] <- TRUE
        uruns <- rle(inverse.rle(uruns))
      }
    }
    ## find events whose length is less than 'mindur' and delete them
    if (mindur > 1) {
      nonev <- with(uruns, values == TRUE & lengths < mindur)
      if (any(nonev)) {
        ## (leave initial period alone)
        nonev[1] <- FALSE
        ## delete short events
        uruns$values[nonev] <- FALSE
        uruns <- rle(inverse.rle(uruns))
      }
    }

    ## TODO: could just return vector index(x) with first item in each group repeated, with NAs for gaps?
    ## - but aggregate.zoo falls over if passed NAs in 'by'.

    ## assign unique numbers to events
    runvals <- uruns$values
    if (all) {
      ## assign numbers to events and inter-events
      ids <- seq_along(runvals)
      uruns$values <- ids
    } else {
      ## assign numbers to events, leave rest as NA
      ids <- seq_len(sum(runvals))
      uruns$values[runvals == TRUE] <- ids
      uruns$values[runvals == FALSE] <- NA
    }
    ## expand back into time series
    ev <- inverse.rle(uruns)
    ## set labels to the index() value of start of each event
    starts <- !duplicated(ev)
    starts[is.na(ev)] <- FALSE
    ## convert to an ordered factor
    attr(ev, "levels") <- as.character(index(x)[starts])
    class(ev) <- c("ordered", "factor")
    ## return events as a zoo with "factor" coredata
    ans <- zoo(ev, index(x))
    if (continue) {
      ans <- na.locf(ans, na.rm = FALSE)
    }
    attr(ans, "thresh") <- thresh
    ans
  }

#' @rdname events
#' @export
eventapply <-
  function(X, events,
           FUN = sum, ..., by.column = TRUE, simplify = TRUE,
           TIMING = c("start", "middle", "end")) {
    FUN <- match.fun(FUN)
    TIMING <- match.arg(TIMING)
    Xnames <- colnames(X)
    if (inherits(events, "zoo") || inherits(events, "ts")) {
      ## merge series together (typically zoo or ts objects)
      cX <- cbind(X, events)
      events <- cX[, ncol(cX)]
      X <- cX[, -ncol(cX)]
      colnames(X) <- Xnames
    }
    events <- coredata(events)
    ## need to handle functions returning vectors as well as scalars
    if (length(dim(X)) == 2) {
      ## multiple series
      if (by.column) {
        ## apply to each column separately
        ans <-
          lapply(as.data.frame(X), function(x) {
            tmp <-
              sapply(split(x, events, drop = TRUE),
                FUN, ...,
                simplify = simplify
              )
            if (is.matrix(tmp)) t(tmp) else tmp
          })
        ans <- as.matrix(do.call("data.frame", ans))
      } else {
        ## pass sub-period of multivariate series to function
        ans <-
          sapply(split(seq_len(NROW(X)), events, drop = TRUE),
            function(ii) FUN(X[ii, , drop = FALSE], ...),
            simplify = simplify
          )
        if (is.matrix(ans)) {
          ans <- t(ans)
        }
      }
    } else {
      ## only one series
      ans <- sapply(split(coredata(X), events, drop = TRUE),
        FUN, ...,
        simplify = simplify
      )
      if (is.matrix(ans)) {
        ans <- t(ans)
      }
    }
    ## reduce to plain vector if only one dimension
    if (NCOL(ans) == 1) {
      ans <- drop(ans)
    }
    ## select time corresponding to each event
    timeIdxFn <- switch(TIMING,
      start = function(x) x[1],
      middle = function(x) x[ceiling(length(x) / 2)],
      end = function(x) x[length(x)]
    )
    ev.index <- unlist(lapply(
      split(seq_len(NROW(X)), events, drop = TRUE),
      timeIdxFn
    ))
    ev.times <- time(X)[ev.index]
    if (simplify && !is.list(ans)) {
      return(zoo(as.matrix(ans), ev.times))
    } else {
      names(ans) <- format(ev.times)
      return(ans)
    }
  }

preInterEventDuration <- function(x) {
  if (inherits(x, "zoo")) {
    stopifnot(is.regular(x, strict = TRUE))
  }
  interCounter <- cumsum(is.na(coredata(x)))
  vals <- tapply(interCounter, coredata(x), FUN = head, 1)
  c(vals[1], diff(vals))
}

#' @rdname events
#' @export
eventinfo <-
  function(X, events,
           FUN = mean, ...) {
    stopifnot(inherits(X, "zoo"))
    xValues <- eventapply(X, events = events, FUN = FUN, ...)
    xLengths <- eventapply(X,
      events = events, FUN = NROW,
      TIMING = "middle", by.column = FALSE
    )
    midTimeComponents <- as.POSIXlt(time(xLengths))
    data.frame(
      Time = time(xValues),
      Month = midTimeComponents$mon + 1,
      Year = midTimeComponents$year + 1900,
      Value = coredata(xValues),
      Duration = coredata(xLengths),
      PreDuration = preInterEventDuration(events)
    )
  }
