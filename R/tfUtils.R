## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


tfParsCheck <-
  function(pars) {
    if (length(pars) == 0) {
      return(TRUE)
    }
    pars <- tfParsConvert(pars, "a,b")
    a <- pars[grep("^a", names(pars))]
    b <- pars[grep("^b", names(pars))]
    stopifnot(length(b) > 0)
    ## check that parameters are real numbers
    stopifnot(all(is.finite(pars)))
    ## check lambda
    if ("lambda" %in% names(pars)) {
      stopifnot(pars["lambda"] >= -1)
      stopifnot(pars["lambda"] <= 0)
    } else {
      ## check that AR component is in the region of stationarity
      ## (copied from stats::arima)
      if (length(a) > 0) {
        if (!all(Mod(polyroot(c(1, -a))) > .95)) {
          stop("AR component not in the region of stationarity")
        }
      }
    }
    ## convert parameters to more intuitive form (for up to 2nd order only)
    #    if (isStandardModelOrder(pars)) {
    #        pars.tv <- abToTauV(pars)
    #        stopifnot(all(is.finite(pars.tv)))
    #        parSymbols <- gsub("_.*", "", names(pars.tv))
    ## check that recessions are real and all components have positive volume
    #        taus <- (parSymbols == "tau")
    #        vs <- (parSymbols == "v")
    #        stopifnot(all(pars.tv[taus] >= 0))
    #        stopifnot(all(pars.tv[vs] >= 0))
    #    }
    return(TRUE)
  }

tfParsConvert <-
  function(pars, form = c("a,b", "tau,v", "alpha,beta")) {
    form <- match.arg(form)
    pars <- unlist(pars)
    ## extract first parts of parameter names (before underscore)
    parSymbols <- gsub("_.*", "", names(pars))

    if (any(c("tau", "v") %in% parSymbols)) {
      if (form == "a,b") {
        return(do.call("tauVToAB", as.list(pars)))
      }
      if (form == "alpha,beta") {
        return(do.call("tauVToAB", c(as.list(pars), list(alphabeta = TRUE))))
      }
      return(pars)
    }
    if (any(c("a", "b") %in% parSymbols)) {
      if (form == "tau,v") {
        return(abToTauV(pars))
      }
      if (form == "alpha,beta") {
        return(abToTauV(pars, alphabeta = TRUE))
      }
      return(pars)
    }
    # warning("unknown parameter set")
    pars
  }

isStandardModelOrder <- function(pars) {
  pars <- unlist(pars)
  ## extract first parts of parameter names (before underscore)
  parSymbols <- gsub("_.*", "", names(pars))
  if (any(c("tau", "v") %in% parSymbols)) {
    return(TRUE)
  }
  ## assuming a,b form now
  a <- pars[parSymbols == "a"]
  b <- pars[parSymbols == "b"]
  a <- stripzeros(a)
  b <- stripzeros(b, up.to = 1)
  n <- length(a)
  m <- length(b) - 1
  ## we can only represent certain a,b forms as tau,v:
  ok <- (n <= 3) && (m <= 3)
  ok <- ok && (m <= n + 1)
  # ok <- ((n <= 2) && (n >= m) && (n - m <= 1))
  return(ok)
}

stripzeros <- function(x, up.to = 0) {
  eps <- sqrt(.Machine$double.eps)
  nonZeros <- which(abs(x) > eps)
  if (length(nonZeros) > 0) {
    up.to <- max(tail(nonZeros, 1), up.to)
  }
  x[seq_len(up.to)]
}

polesToAr <- function(poles) {
  ar <- numeric()
  if (length(poles) > 0) {
    ar[1] <- sum(poles)
    if (length(poles) == 2) {
      ar[2] <- -prod(poles)
    }
    if (length(poles) == 3) {
      ar[2] <- -(prod(poles[-3]) + prod(poles[-2]) + prod(poles[-1]))
      ar[3] <- prod(poles)
    }
  }
  ar
}

arToPoles <- function(ar) {
  aa <- c(rev(-ar), 1)
  aastab <- polystab(aa)
  poles <- polyroot(aastab)
  eps <- sqrt(.Machine$double.eps)
  if (!any(abs(Im(poles)) > eps)) {
    poles <- Re(poles)
  }
  poles
}

abToTauV <-
  function(theta, a, b, alphabeta = FALSE) {
    if (missing(theta)) {
      theta <- c(a, b)
      n <- length(a)
      m <- length(b) - 1
      a_names <- paste("a", 1:n, sep = "_")
      b_names <- paste("b", 0:m, sep = "_")
      names(theta) <- c(a_names, b_names)
    } else {
      if (!missing(a) || !missing(b)) {
        stop("give either theta or a,b")
      }
    }
    if (!isStandardModelOrder(theta)) {
      stop("can not represent this model order as tau,v")
    }
    ## extract first parts of parameter names (before underscore)
    pars <- unlist(theta)
    parSymbols <- gsub("_.*", "", names(pars))
    a <- pars[parSymbols == "a"]
    b <- pars[parSymbols == "b"]
    ## keep any extra parameters
    extraPars <- pars[(parSymbols %in% c("a", "b")) == FALSE]
    a <- unname(a)
    b <- unname(b)

    stopifnot(length(a) <= 3)
    length(a) <- 3
    if (length(b) < 3) length(b) <- 3
    a[!is.finite(a)] <- 0
    b[!is.finite(b)] <- 0

    eps <- sqrt(.Machine$double.eps)

    aa <- c(rev(-a), 1)
    aastab <- polystab(aa)
    poles <- polyroot(aastab)
    # poles <- polyroot(c(rev(-a), 1))
    if (any(Re(poles) < -eps)) {
      warning("negative poles detected")
      return(NA)
    }
    if (any(abs(Im(poles)) > eps)) {
      warning("complex poles detected")
    }
    alpha_s <- max(Re(poles))
    alpha_q <- median(Re(poles))
    alpha_3 <- min(Re(poles))
    # alpha_s <- max(alpha_s, 0)
    # alpha_q <- max(alpha_q, 0)
    # alpha_3 <- max(alpha_3, 0)

    series <- NA
    if ("series" %in% names(pars)) {
      series <- round(pars[["series"]])
    }

    if (is.na(series)) {
      series <- 0

      if (alpha_q == 0) {
        ## first-order model
        series <- 0
      } else if (alpha_3 == 0) {
        ## second-order model
        if (all(b[-1] == 0)) {
          series <- 1
        }
        if (abs(alpha_s - alpha_q) < eps) {
          series <- 1
        }
      } else {
        ## third-order model

        ## TODO: how to distinguish series=1 from series=0?
        # if (FALSE) {
        alpha <- c(alpha_s, alpha_q, alpha_3)
        dup <- (duplicated(signif(alpha, 5)) |
          duplicated(signif(alpha, 5), fromLast = TRUE))
        if (any(dup)) {
          series <- 1
          ## by convention, components q & 3 are in series
          adup <- head(alpha[dup], 1)
          alpha_q <- alpha_3 <- adup
          alpha_s <- head(alpha[!dup], 1)
        }
        # }
        if (all(b[-(1:2)] == 0)) {
          series <- 2
        }
        if (all(b[-1] == 0)) {
          series <- 3
        }
      }
    }

    ## TODO: use residues function from Octave?

    if (series == 0) {
      ## all components in parallel
      if (alpha_s == 0) {
        beta_s <- b[1]
        beta_q <- 0
        beta_3 <- 0
      } else if ((alpha_q == 0) && (alpha_3 == 0)) {
        ## first-order model
        beta_q <- -b[2] / alpha_s
        beta_s <- b[1] - beta_q
        beta_3 <- 0
      } else {
        ## second or third order
        AA <- rbind(
          c(1, 1, 1),
          -c(alpha_q + alpha_3, alpha_s + alpha_3, alpha_s + alpha_q),
          c(alpha_q * alpha_3, alpha_s * alpha_3, alpha_s * alpha_q)
        )
        beta <- solve(AA, b[1:3])
        beta_s <- beta[1]
        beta_q <- beta[2]
        beta_3 <- beta[3]
      }
    } else {
      if (alpha_3 == 0) {
        ## second-order model in series
        ## take gain of slow component to be 1
        beta_s <- (1 - alpha_s)
        # beta_q <- (1 - alpha_q)
        beta_q <- b[1] / beta_s
        ## TODO: what if b[2] != 0? (equal roots case)
        beta_3 <- 0
      } else {
        ## third-order model
        if (series == 1) {
          ## one component in parallel with two in series
          ## (q & 3 are in series; s in parallel)
          beta_s <- b[3] / (alpha_q * alpha_3)
          beta_q <- sqrt((-b[2] - beta_s * (alpha_q + alpha_3)) /
            alpha_s * (b[1] - beta_s))
          beta_3 <- (b[1] - beta_s) / beta_q
        } else if (series == 2) {
          ## one component in series with two in parallel
          ## (3 in series; s & q in parallel)
          ## apply gain to v_3 only:
          gain <- sum(b) / (1 - sum(a))
          beta_3 <- (1 - alpha_3) * gain
          beta_q <- ((b[2] + alpha_q * b[1]) /
            (beta_3 * (alpha_q - alpha_s)))
          beta_s <- (b[1] / beta_3) - beta_q
        } else if (series == 3) {
          ## three components in series
          beta_s <- (1 - alpha_s)
          beta_q <- (1 - alpha_q)
          beta_3 <- (1 - alpha_3)
          ## apply gain to v_3 only:
          beta_3 <- beta_3 * (b[1] / (beta_s * beta_q * beta_3))
        } else {
          stop("unrecognised values of 'series': ", series)
        }
      }
    }

    if (any(!is.finite(c(beta_s, beta_q, beta_3)))) {
      return(NA)
    }

    if (alphabeta) {
      return(c(
        if (alpha_s != 0) c(alpha_s = alpha_s),
        if (alpha_q != 0) c(alpha_q = alpha_q),
        if (alpha_3 != 0) c(alpha_3 = alpha_3),
        if (TRUE) c(beta_s = beta_s),
        if (beta_q != 0) c(beta_q = beta_q),
        if (beta_3 != 0) c(beta_3 = beta_3),
        extraPars,
        if (series > 0) c(series = series)
      ))
    }

    ## convert from 'alpha' and 'beta' to 'tau' and 'v'
    tau_q <- -1 / log(alpha_q)
    tau_s <- -1 / log(alpha_s)
    tau_3 <- -1 / log(alpha_3)
    v_s <- beta_s / (1 - alpha_s)
    v_q <- beta_q / (1 - alpha_q)
    v_3 <- beta_3 / (1 - alpha_3)

    with_tq <- ((tau_q != 0) || (("lambda" %in% parSymbols)))
    with_vq <- ((v_q != 0) && (("lambda" %in% parSymbols) == FALSE))
    c(
      if (tau_s != 0) c(tau_s = tau_s),
      if (with_tq) c(tau_q = tau_q),
      if (tau_3 != 0) c(tau_3 = tau_3),
      if (TRUE) c(v_s = v_s),
      if (with_vq) c(v_q = v_q),
      if (v_3 != 0) c(v_3 = v_3),
      extraPars,
      if (series > 0) c(series = series)
    )
  }

tauVToAB <-
  function(tau_s = 0, tau_q = 0, tau_3 = 0,
           v_s = 1, v_q = NA, v_3 = 0,
           ..., series = 0, alphabeta = FALSE) {
    series <- round(series)
    stopifnot(series %in% 0:3)
    if ((series > 1) && (v_3 == 0)) {
      stop("series > 1 is only valid with v_3 > 0")
    }
    ## fix ridiculous values (avoid numerical problems when converting back)
    tau_s <- min(tau_s, 1e10)
    tau_q <- min(tau_q, 1e4)
    tau_3 <- min(tau_3, 1e4)
    if (tau_s < (-1 / log(1e-6))) tau_s <- 0
    if (tau_q < (-1 / log(1e-6))) tau_q <- 0
    if (tau_3 < (-1 / log(1e-6))) tau_3 <- 0
    ## default value of v_q depends on 'series'
    ## (by convention only, this defines configuration of the system)
    if (is.na(v_q)) {
      if (series == 0) {
        ## all in parallel
        v_q <- 1 - v_s - v_3
      }
      if (series == 1) {
        v_q <- 1
      }
      if (series == 2) {
        v_q <- 1 - v_s
      }
      if (series == 3) {
        v_q <- 1
      }
      v_q <- max(0, min(1, v_q))
    }
    ## convert from 'tau' and 'v' to 'alpha' and 'beta'
    alpha_s <- exp(-1 / tau_s)
    alpha_q <- exp(-1 / tau_q)
    alpha_3 <- exp(-1 / tau_3)
    beta_s <- v_s * (1 - alpha_s)
    beta_q <- v_q * (1 - alpha_q)
    beta_3 <- v_3 * (1 - alpha_3)
    if (alpha_s == 1) {
      beta_s <- v_s
    }
    if (alphabeta) {
      return(c(
        alpha_s = alpha_s,
        alpha_q = alpha_q,
        alpha_3 = alpha_3,
        beta_s = beta_s,
        beta_q = beta_q,
        beta_3 = beta_3,
        ...
        # if (series > 0) c(series = series)
      ))
    }

    ## convert from 'alpha' and 'beta' to 'a' and 'b'
    a <- b <- numeric(3)
    a[1] <- alpha_s + alpha_q + alpha_3
    a[2] <- -(alpha_s * alpha_q +
      alpha_s * alpha_3 +
      alpha_q * alpha_3)
    a[3] <- alpha_s * alpha_q * alpha_3

    if (series == 0) {
      ## all stores in parallel
      b[1] <- beta_s + beta_q + beta_3
      b[2] <- -(beta_s * (alpha_q + alpha_3) +
        beta_q * (alpha_s + alpha_3) +
        beta_3 * (alpha_s + alpha_q))
      b[3] <- (beta_s * alpha_q * alpha_3 +
        beta_q * alpha_s * alpha_3 +
        beta_3 * alpha_s * alpha_q)
    } else {
      if (beta_3 == 0) {
        ## second-order model in series
        b[1] <- beta_s * beta_q
        b[-1] <- 0
      } else {
        ## third-order model
        if (series == 1) {
          ## two components in series and one in parallel
          ## (q & 3 are in series; s in parallel)
          b[1] <- beta_q * beta_3 + beta_s
          b[2] <- -(beta_q * beta_3 * alpha_s +
            beta_s * (alpha_q + alpha_3))
          b[3] <- beta_s * alpha_q * alpha_3
        } else if (series == 2) {
          ## one component in series with two in parallel
          ## (3 in series; s & q in parallel)
          b[1] <- (beta_s + beta_q) * beta_3
          b[2] <- -(beta_s * alpha_q + beta_q * alpha_s) * beta_3
          b[-(1:2)] <- 0
        } else if (series == 3) {
          ## three components in series
          b[1] <- beta_s * beta_q * beta_3
          b[-1] <- 0
        } else {
          stop("unrecognised values of 'series': ", series)
        }
      }
    }

    a <- stripzeros(a)
    b <- stripzeros(b, up.to = 1)
    n <- length(a)
    m <- length(b) - 1
    theta <- c(a, b, ...) # , if (series > 0) c(series = series))
    if (n > 0) {
      names(theta)[1:n] <- paste("a", 1:n, sep = "_")
    }
    names(theta)[(n + 1) + 0:m] <- paste("b", 0:m, sep = "_")
    theta
  }


stabiliseAR <- function(ar) {
  if (length(ar) < 1) {
    return(ar)
  }
  a <- c(rev(-ar), 1)
  aa <- polystab(a)
  ## desired form has a[n] = 1
  aa <- rev(aa)
  aa[-1] <- aa[-1] / aa[1]
  arnew <- -aa[-1]
  return(arnew)
}

## polynomial roots stabilisation
## based on code from Peter Young 2009-11-06
polystab <- function(a) {
  if (length(a) <= 1) {
    return(a)
  }
  v <- polyroot(a)
  i <- which(v != 0)
  vs <- 0.5 * (sign(abs(v[i]) - 1) + 1)
  v[i] <- (1 - vs) * v[i] + vs / Conj(v[i])
  b <- a[a != 0]
  suppressWarnings({
    b <- b[1] * coef(poly.calc(v))
  })
  if (!any(Im(a) != 0)) {
    b <- Re(b)
  }
  b
}


ssg.tf.coef <-
  function(theta, ...) {
    if (length(theta) == 0) {
      return(NA_real_)
    }
    ## loss term interferes with gain calculation
    if ("loss" %in% names(theta)) {
      if (theta[["loss"]] != 0) {
        return(NA_real_)
      }
    }
    ## in lambda model, gain is defined to be 1
    ## (simulation function constrains total volume to 1)
    if ("lambda" %in% names(theta)) {
      if (theta[["lambda"]] != 0) {
        return(1)
      }
    }
    a <- theta[grep("^a", names(theta))]
    b <- theta[grep("^b", names(theta))]
    if (length(b) == 0) {
      return(NA_real_)
    }
    sum(b) / (1 - sum(a))
  }

normalise.tf.coef <-
  function(theta, ...) {
    gain <- ssg.tf.coef(theta, ...)
    if (is.na(gain) || (gain == 1) || (gain == 0)) {
      return(theta)
    }
    ib <- grep("^b", names(theta))
    b <- theta[ib]
    theta[ib] <- b / gain
    theta
  }

describeTF <- function(theta, ...) {
  pars <- theta
  if (length(pars) == 0) {
    return(NA)
  }
  pars <- tfParsConvert(pars, "a,b")
  ## model order
  parSymbols <- gsub("_.*", "", names(pars))
  a <- pars[parSymbols == "a"]
  b <- pars[parSymbols == "b"]
  a <- stripzeros(a)
  b <- stripzeros(b, up.to = 1)
  n <- length(a)
  m <- length(b) - 1
  if ((n == 0) && (m < 0)) {
    return(NA)
  }
  ans <- NULL
  if (n == 0) {
    ans <- "instantaneous"
  } else if (n == 1) {
    if (m == 0) ans <- "single store"
    if (m == 1) ans <- "single store + instantaneous in parallel"
    if (m >= 2) ans <- "single store + complex MA component"
  } else if (n == 2) {
    if (m == 0) ans <- "S * Q (two stores in series)"
    if (m == 1) ans <- "S + Q (two stores in parallel)"
    if (m == 2) ans <- "S + Q + inst. (three in parallel)"
    if (m >= 3) ans <- "S + Q + complex MA component"
  } else if (n == 3) {
    series <- 0
    if ("series" %in% names(theta)) {
      series <- round(theta[["series"]])
    } else {
      partv <- try(tfParsConvert(pars, "tau,v"))
      if ("series" %in% names(partv)) {
        series <- round(partv[["series"]])
      }
    }
    if (series == 0) ans <- "S + Q + 3 (three in parallel)"
    if (series == 1) ans <- "(Q * 3) + S (two in series, one in parallel)"
    if (series == 2) ans <- "(S + Q) * 3 (two in parallel, one in series)"
    if (series == 3) ans <- "S * Q * 3 in series"
    if ((m == 0) != (series == 3)) {
      warning("m == ", m, " but series == ", series)
    }
  } else {
    ans <- "complex (high order?)"
  }
  if (length(a > 1)) {
    poles <- arToPoles(a)
    ans <- paste(ans, "\n", "    Poles:", toString(signif(poles, 4)), sep = "")
  }
  ans
}
