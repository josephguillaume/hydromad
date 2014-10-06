##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

findThresh <-
    function(x, thresh = NA, ## ignored
             n, within = (n %/% 20) + 1,
             mingap = 1, mindur = 1, 
             below = FALSE, ...,
             trace = FALSE, optimize.tol = 0.1)
{
    stopifnot(n > 0)
    n <- round(n)
    stopifnot(n * mindur + n * mingap <= NROW(x))
    x <- coredata(x)
    if (below) {
        x <- -x
    }
    ## return (difference from 'n' of) number of events for 'thresh'
    nDiffForThresh <- function(thresh) {
        ev <- eventseq(x, thresh = thresh, mingap = mingap,
                       mindur = mindur, ...)
        newn <- nlevels(ev)
        if (trace)
            message(sprintf("thresh = %.3f, n = %d", thresh, newn))
        abs(newn - n)
    }

    ## do a quick run with a rough guess for the range
    rng <- quantile(x, (1 - (n * mindur / NROW(x))) ^ c(1, 3),
                     na.rm = TRUE, names = FALSE)
    if (trace)
        message("initial range guess: ", toString(signif(rng, 4)))
    res <- optimize(nDiffForThresh, interval = rng,
                    tol = optimize.tol)
    if (res$objective > within) {
        if (trace)
            message("initial run off by ", res$objective)
        ## expand possible range; only exclude minimum and maximum values
        rng2 <- quantile(x, range(ppoints(length(x))),
                        na.rm = TRUE, names = FALSE)
        res2 <- optimize(nDiffForThresh, interval = rng2,
                         tol = optimize.tol)
        if (res2$objective < res$objective)
            res <- res2
    }
    if (res$objective > within)
        warning("actual number of events differs from target (", n, ") by ",
                res$objective)
    thresh <- res$minimum
    if (below)
        thresh <- -thresh
    return(thresh)
}

eventseq <-
    function(x, thresh = 0, mingap = 1, mindur = 1, extend = 0,
             inthresh = thresh, inx = x, indur = 1,
             below = FALSE, all = FALSE, continue = FALSE,
             n = NA)
{
    if (!is.na(n) && (n > 0)) {
        if (thresh > 0)
            warning("'thresh' argument ignored since 'n' given")
        ccall <- match.call()
        ccall[[1]] <- quote(findThresh)
        thresh <- eval.parent(ccall)
    }
    ## check for simple case:
    inthreshGiven <-
        ( ((!missing(inthresh) || !missing(inx)) &&
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
        if (length(thresh) > 1)
            threshmat <- matrix(thresh, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
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
            if (length(inthresh) > 1)
                inthreshmat <- matrix(inthresh, nrow = nrow(inx), ncol = ncol(inx), byrow = TRUE)
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
    if (continue)
        ans <- na.locf(ans, na.rm = FALSE)
    attr(ans, "thresh") <- thresh
    ans
}

eventapply <-
    function(X, events,
             FUN = sum, ..., by.column = TRUE, simplify = TRUE,
             TIMING = c("start", "middle", "end"))
{
    FUN <- match.fun(FUN)
    TIMING <- match.arg(TIMING)
    Xnames <- colnames(X)
    if (inherits(events, "zoo") || inherits(events, "ts")) {
        ## merge series together (typically zoo or ts objects)
        cX <- cbind(X, events)
        events <- cX[,ncol(cX)]
        X <- cX[,-ncol(cX)]
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
                               FUN, ..., simplify = simplify)
                    if (is.matrix(tmp)) t(tmp) else tmp
                })
            ans <- as.matrix(do.call("data.frame", ans))
        } else {
            ## pass sub-period of multivariate series to function
            ans <-
                sapply(split(seq_len(NROW(X)), events, drop = TRUE),
                       function(ii) FUN(X[ii,,drop=FALSE], ...),
                       simplify = simplify)
            if (is.matrix(ans))
                ans <- t(ans)
        }

    } else {
        ## only one series
        ans <- sapply(split(coredata(X), events, drop = TRUE),
                      FUN, ..., simplify = simplify)
        if (is.matrix(ans))
            ans <- t(ans)
    }
    ## reduce to plain vector if only one dimension
    if (NCOL(ans) == 1)
        ans <- drop(ans)
    ## select time corresponding to each event
    timeIdxFn <- switch(TIMING,
                        start = function(x) x[1],
                        middle = function(x) x[ceiling(length(x)/2)],
                        end = function(x) x[length(x)])
    ev.index <- unlist(lapply(split(seq_len(NROW(X)), events, drop = TRUE),
                       timeIdxFn))
    ev.times <- time(X)[ev.index]
    if (simplify && !is.list(ans)) {
        return(zoo(as.matrix(ans), ev.times))
    } else {
        names(ans) <- format(ev.times)
        return(ans)
    }
}

preInterEventDuration <- function(x)
{
    if (inherits(x, "zoo"))
        stopifnot(is.regular(x, strict = TRUE))
    interCounter <- cumsum(is.na(coredata(x)))
    vals <- tapply(interCounter, coredata(x), FUN = head, 1)
    c(vals[1], diff(vals))
}

eventinfo <-
    function(X, events,
             FUN = mean, ...)
{
    stopifnot(inherits(X, "zoo"))
    xValues <- eventapply(X, events = events, FUN = FUN, ...)
    xLengths <- eventapply(X, events = events, FUN = NROW,
                           TIMING = "middle", by.column = FALSE)
    midTimeComponents <- as.POSIXlt(time(xLengths))
    data.frame(Time = time(xValues),
               Month = midTimeComponents$mon + 1,
               Year = midTimeComponents$year + 1900,
               Value = coredata(xValues),
               Duration = coredata(xLengths),
               PreDuration = preInterEventDuration(events))
}
