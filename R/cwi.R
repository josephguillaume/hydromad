## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## IHACRES Catchment Wetness Index (CWI) model.
cwi.sim <-
  function(DATA,
           tw, f = 0, scale,
           l = 0, p = 1,
           t_ref = hydromad.getOption("cwi")$t_ref,
           s_0 = 0,
           return_state = FALSE) {
    if (NCOL(DATA) > 1) stopifnot("P" %in% colnames(DATA))
    if (f > 0) stopifnot("E" %in% colnames(DATA))
    ## special value scale = NA used for initial run for scaling
    if (is.na(scale)) {
      scale <- 1
    }
    ## check values
    stopifnot(tw >= 0)
    stopifnot(f >= 0)
    stopifnot(scale >= 0)
    stopifnot(l >= 0)
    stopifnot(p > 0)
    stopifnot(length(t_ref) == 1)

    inAttr <- attributes(if (NCOL(DATA) > 1) DATA[, 1] else DATA)
    DATA <- as.ts(DATA)
    P <- if (NCOL(DATA) > 1) DATA[, "P"] else DATA

    ## compute soil moisture index s
    if (f == 0) {
      ## this is an invariant drying rate model
      w <- 1 - 1 / max(tw, 1)
      w <- rep(w, length(P))
      ## use filter_tv rather than filter because it maintains
      ## state over NA gaps (filter resets to 0)
      # s <- filter(P, w, method="recursive", init = s_0)
    } else {
      tw_k <- tw * exp(-0.062 * f * (DATA[, "E"] - t_ref))
      tw_k <- pmax(tw_k, 1)
      w <- 1 - 1 / tw_k
    }
    s <- filter_tv(P, w, init = s_0)
    ## compute effective rainfall U
    U <- (pmax(s - l, 0)^p) * P
    ## TODO: test mass balance scaling with p > 1
    U <- (scale^p) * U

    attributes(U) <- inAttr
    ans <- U
    if (return_state) {
      attributes(s) <- attributes(w) <- attributes(U)
      ans <- cbind(U = U, s = s, w = w)
    }
    return(ans)
  }

cwi.ranges <- function() {
  list(
    tw = c(0, 100),
    f = c(0, 8),
    scale = NA_real_,
    l = 0,
    p = 1,
    t_ref = 20
  )
}

absorbScale.hydromad.cwi <- function(object, gain, ...) {
  parname <- "scale"
  if (gain <= 0) {
    return(NULL)
  }
  coeff <- coef(object, which = "sma")
  scale <- coeff[[parname]]
  ## we only want to do this when scale is NA (special value)
  if (is.null(scale) || !is.na(scale)) {
    return(NULL)
  }
  scale <- 1
  p <- coeff[["p"]]
  scale <- scale * (gain^p)
  scale <- max(scale, 0)
  # object <- update(object, scale = scale, and.rescale = FALSE)
  ## may be significatly faster:
  object$parlist[[parname]] <- scale
  object$call[[parname]] <- signif(scale, 6)
  object$U <- (scale^p) * object$U
  object
}
