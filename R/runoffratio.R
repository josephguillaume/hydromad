## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

runoffratio.sim <-
  function(DATA, width = 30, kernel = 2, sides = 2, rrthresh = 0,
           qlag = 0, scale = 1, return_state = FALSE) {
    stopifnot(c("P", "Q") %in% colnames(DATA))
    P <- DATA[, "P"]
    Q <- DATA[, "Q"]
    ## special value scale = NA used for initial run for scaling
    if (is.na(scale)) {
      scale <- 1
    }
    ## check values
    stopifnot(kernel >= 1)
    ## synchronise Q lagged by qlag
    Q <- shiftWindow(Q, round(qlag), and.lag = TRUE)
    ## compute effective rainfall U
    ## estimate U as scaled P, scaled in a moving window
    ## (runoff coefficient)
    width <- round(width)
    sm <- simpleSmoothTs(cbind(Q, P), width = width, c = kernel, sides = sides)
    rr <- sm[, "Q"] / sm[, "P"]
    rr <- na.locf(na.locf(rr, na.rm = FALSE), fromLast = TRUE, na.rm = FALSE)
    rr[!is.finite(rr)] <- 0
    rr[rr < rrthresh] <- 0
    U <- scale * P * rr

    ans <- U
    if (return_state) {
      ans <- cbind(U = U, rr = rr)
    }
    return(ans)
  }

runoffratio.ranges <- function() {
  list(
    rrthresh = c(0, 0.2),
    scale = NA
  )
}

absorbScale.hydromad.runoffratio <-
  function(object, gain, ...) {
    absorbScale.hydromad.scalar(object, gain, parname = "scale")
  }
