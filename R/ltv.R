
# Interactive SRW estimation with change detection
# based on version 15.3.07 of John Norton's MatLab code
# R port Copyright (c) 2007 Felix Andrews <felix@nfrac.org>

# Function  arguments:
# y: output vector
# u: input vector or matrix with observations corresponding to y
# na: no. of denominator terms
# nb: no. of numerator terms, i.e. order+1 -- for each column of u (will be repeated if necessary)
# nc: noise-model order
# delay: dead times; first non-zero unit-pulse response is at time delay+1 -- for each column of u (will be repeated if necessary)
# qa: variance of random walk of a parameters (repeated if necessary)
# qb: variance of random walk of b parameters (repeated if necessary)
# qd: variance of random walk of offset
# Q: covariance matrix of random walk of all parameters
#    by default this is diag(c(rep(q_a, length=na), rep(q_b, length=sum(nb)), q_d))

# need a scheme for iterating and looking at diagnostics
# to choose appropriate model order and degree of parameter time variation
# (Norton and Chanat paper: reduce mse but keep v_e_ratio very little above 1
#                           AND assess credibility of parameter variation
# "when na is too high, one or more poles will become unrealistic (e.g. complex)
#  or will have a very small associated modal amplitude."
# "the noise order nc is a compromise between flexibility and parsimony."


# The only thing I have not implemented is calculation of residues. Poles are easy enough with the "polyroot" function, but R does not seem to have an equivalent to the Matlab "residue" function. Since I have no idea what "residues" are, I can not guess how to implement it.

# U <- ts(rpois(200, lambda=0.25))
# Q <- filter(0.5 * U, filter=0.5, method="recursive")
# Q <- Q + rnorm(200, sd=0.1)
# Q <- pmax(Q, 0)
# xyplot(ts.union(Q, U))
# foo <- ltv(Q, U, na=1, nb=1, nc=1)
# xyplot(foo$filter)
# xyplot(ts.union(foo$y, foo$sim.y, foo$residuals, foo$innov))
# TODO: multi-variate example etc

# TODO: doc!


# library(dlm)
#
# QXX <- ts.intersect(Q, lag(Q, -1), lag(Q, -2), lag(U, -1), lag(U, -2))
# mod <- dlmModReg(QXX[,-1], addInt=FALSE, dW=c(0,0,1e-5,1e-5))
# f <- dlmFilter(QXX[,1], mod)
# xyplot(window(f$m, start=start(U)[1]+10))
# s <- dropFirst(dlmSmooth(f)$s)
# xyplot(s)
# sim-mode output... (use s and U)
# ysim <- ts(0, start=start(s), end=end(s))
# for (t in seq(3,NROW(s))) {
# 	x_t <- c(ysim[t-1], ysim[t-2], U[t-1], U[t-2])
# 	ysim[t] <- s[t,] %*% x_t
# }
# resid <- Q - ysim
# xyplot(ts.union(y=f$y, sim=ysim,
# 	sim_residuals=resid, forecast_residuals=residuals(f,sd=F)))
# xyplot(ts.union(y=f$y, sim=ysim, resid=resid, forecast_err=residuals(f,sd=F), s)

### TODO: does any indexing etc assume that time series start at 1?

# alternative interface: function(formula, data, na.action=na.fail, ...
ltv <- function(y, u, na = 3, nb = 3, nc = 3, delay = 0, q_a = 0, q_b = 0, q_c = 0, q_d = 0, Q = NULL) {
  stopifnot(is.numeric(y))
  stopifnot(is.numeric(u) || is.data.frame(u))
  stopifnot(NROW(u) == NROW(y))
  # TODO: what about missing values?
  # mf <- model.frame(formula, data=data, na.action=na.action)
  # mt <- attr(mf, "terms")
  # y <- as.ts(model.response(mf, "numeric"))
  # x <- as.ts(model.matrix(mt, mf))
  if (!is.ts(y)) y <- as.ts(y)
  if (!is.matrix(u)) {
    tmp <- tsp(u)
    u <- as.matrix(u)
    tsp(u) <- tmp
    colnames(u) <- "u"
  }
  if (!is.ts(u)) u <- as.ts(u)
  nt <- nrow(u)
  nu <- ncol(u)
  nb <- rep(nb, length = nu)
  nab <- na + sum(nb)
  np <- nab + nc + 1
  # indices to select subsets of parameters
  a.index <- 1:na
  b.index <- (na + 1):(nab)
  abd.index <- c(1:nab, np) # a's, b's and d but not c's
  # parameter names
  xnames <- rep(c("a", "b", "c"), times = c(na, sum(nb), nc))
  xnames <- paste(xnames, sequence(c(na, nb, nc)), sep = "")
  if (nu > 1) {
    xnames[b.index] <- paste(xnames[b.index],
      rep(colnames(u), each = nb),
      sep = "_"
    )
  }
  xnames <- c(xnames, "d")
  # q is diag(Q), gives the random walk variance for each parameter
  q_a <- rep(q_a, length = na)
  q_b <- rep(q_b, length = sum(nb))
  q_c <- rep(q_a, length = nc)
  q <- c(q_a, q_b, q_c, q_d)
  if (is.null(Q)) Q <- diag(q)
  dimnames(Q) <- list(xnames, xnames)
  P <- diag(1e8, np) # normalised P
  dimnames(P) <- list(xnames, xnames)
  hasOffset <- (Q["d", "d"] > 0)
  if (!hasOffset) P["d", "d"] <- 0
  # t0: allows for dead time down to 0 (instant response start)
  delay <- rep(delay, length = nu)
  t0 <- max(c(na, nb + delay, nc)) + 2
  # construct regressor matrix
  H <- list(y = ts.lagbind(y, lags = -(1:na)))
  H[colnames(u)] <- lapply(1:nu, function(i) {
    ts.lagbind(u[, i], lags = -delay[i] - (0:(nb[i] - 1)))
  })
  H <- do.call(ts.union, H)
  H <- window(H, start = start(y), end = end(y), extend = TRUE) # 1:nt ?
  e <- ts(0, start = start(y), end = end(y))
  v <- ts(0, start = start(y), end = end(y))

  ## Construct filtered estimates...

  # time series of filtered parameter estimates
  x.f <- ts(matrix(0, ncol = np), start = start(y), end = end(y))
  colnames(x.f) <- xnames
  PS <- ts(matrix(0, ncol = np), start = start(y), end = end(y))
  colnames(PS) <- xnames
  # initial d estimate = o/p. as a, b, c zero
  if (Q["d", "d"] != 0) {
    x.f[1:(t0 - 1), "d"] <- y[1:(t0 - 1)]
  }
  dn <- 0
  message("\n", "Filtering... ")
  # TODO: write this loop in C
  for (t in t0:nt) {
    # first update is from t0-1 to t0, last to nt
    h <- c(H[t, ], e[t - seq(len = nc)], 1) # fresh regressors, not updated ones
    h[1:na] <- h[1:na] - x.f[t - (1:na), "d"] # completes h[t]
    P <- P + Q # time update (time update of theta does nothing)
    g <- drop(P %*% h) # g(t) of change-detection notes
    dn <- drop(1 + g %*% h)
    v[t] <- y[t] - x.f[t - 1, ] %*% h # innovation t|t-1
    e[t] <- v[t] / dn # residual
    # observation update
    x.f[t, ] <- x.f[t - 1, ] + g * v[t] / dn
    P <- P - g %o% g / dn
    PS[t, ] <- g
  }
  message("Filtering finished.")

  # (sigma^2)(capital lambda) of change-detection notes
  # final condition for sigma^2*LAM recursion
  N <- h %o% h / dn

  ## Construct smoothed estimates...

  # time series of smoothed parameter estimates
  x.s <- ts(matrix(0, ncol = np), start = start(y), end = end(y))
  colnames(x.s) <- xnames
  x.s[nt, ] <- x.f[nt, ]
  dstat <- ts(0, start = start(y), end = end(y)) # change detection statistic
  delta <- ts(matrix(0, ncol = np), start = start(y), end = end(y)) # opt. estimate of state change
  colnames(delta) <- xnames
  cno <- ts(0, start = start(y), end = end(y)) # condition number of N(t)
  mu <- rep(0, np)

  message("\n", "Smoothing... ")
  # TODO: write this loop in C
  for (t in nt:t0) {
    h <- c(H[t, ], e[t - seq(len = nc)], 1) # fresh regressor vector
    h[1:na] <- h[1:na] - x.f[t - (1:na), "d"] # completes h[t]
    gcd <- PS[t, ] # g(t) of change-detection notes (column)
    dn <- drop(1 + gcd %*% h) # denominator for time t

    if (t < nt) {
      fcd <- drop(N %*% gcd) # f(t) of change-detection notes
      N <- (h %o% h) / dn + N -
        (h %o% fcd + fcd %o% h) / dn +
        drop(gcd %*% fcd) * h %o% h / (dn^2) # N(t)
      cno[t] <- kappa(N)
      if (cno[t] < 1e14) {
        delta[t, ] <- solve(-N, mu) # opt. estimate of state change at t
        dstat[t] <- -mu %*% delta[t, ] # mu(t) here; dstat is d*sigma^2
      } else {
        delta[t, ] <- NA
        dstat[t] <- NA
      }
    }
    w <- drop(v[t] + gcd %*% mu)
    mu <- mu - w * h / dn # mu(t) to mu(t-1)
    x.s[t - 1, ] <- x.s[t, ] + Q %*% mu # could use diag(Q)
    v[t] <- y[t] - x.s[t - 1, ] %*% h # innovation
    e[t] <- y[t] - x.s[t, ] %*% h # residual
  }
  message("Smoothing finished.")

  # generate simulation-mode o/p. errors
  ym <- ts(0, start = start(y), end = end(y))
  # initialise regressor vector at time t0-1
  # noise model is omitted
  hs <- y[t0 - 1 - 1:na] - x.s[t0 - 1, "d"] # output - offset
  hs[b.index] <- H[t0 - 1, b.index]
  hs <- c(hs, 1)
  ym[1:(t0 - 1)] <- y[1:(t0 - 1)]
  # TODO: write this loop in C
  for (t in t0:nt) {
    xs <- x.s[t, abd.index] # smoothed a's, b's, d
    if (na > 1) hs[2:na] <- hs[1:(na - 1)]
    hs[1] <- ym[t - 1] - x.s[t - 1, "d"] # model (smoothed) y, not observed y
    hs[b.index] <- H[t, b.index]
    ym[t] <- hs %*% xs
  }

  mean.coef <- colMeans(x.s[, abd.index])
  sd.coef <- apply(x.s[, abd.index], COLS <- 2, sd) # , na.rm=TRUE)
  # apply(x.s[, abd.index], COLS<-2, median, na.rm=TRUE)

  obj <- list(
    call = match.call(),
    na = na, nb = nb, nc = nc, delay = delay, t0 = t0,
    Q = Q, hasOffset = hasOffset,
    filter = x.f[, abd.index],
    smooth = x.s[, abd.index],
    coefficients = mean.coef,
    sd.coefficients = sd.coef,
    # filter=list(
    # 	a=x.f[, a.index],
    # 	b=x.f[, b.index],
    # 	d=x.f[, "d"]
    # ),
    # smooth=list(
    # 	a=x.s[, a.index],
    # 	b=x.s[, b.index],
    # 	d=x.s[, "d"]
    # ),
    y = y,
    # sim.y=ym,
    fitted.values = ym,
    ### residuals=e, # TODO: which residuals are appropriate here??
    residuals = y - ym,
    error = e,
    innovations = v,
    sigma2 = var(v), # ? - for tsdiag
    P = P,
    dstat = dstat,
    delta = delta,
    cno = cno
  )
  class(obj) <- "ltv"
  obj
}

ts.lagbind <- function(x, lags, all = TRUE, dframe = FALSE) {
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

ssgain.ltv <- function(object, ...) {
  if (any(object$smooth$d != 0)) warning("offset is non-zero, so gain will be distorted")
  ans <- rowSums(object$smooth$b) / (1 - rowSums(object$smooth$a))
  mostattributes(ans) <- attributes(object$smooth$a)
  ans
}

poles.ltv <- function(object, step = 10) {
  with(object, {
    steps <- seq(t0, nrow(smooth), by = step)
    nstep <- length(steps)
    A <- cbind(1, -smooth[steps, 1:na])
    # TODO: handle multiple inputs
    B <- smooth[steps, na + 1:nb]
    poles <- ts(matrix(NA, ncol = na, nrow = nstep),
      start = start(smooth), deltat = deltat(smooth) * step
    )
    poles <- apply(A, ROWS <- 1, polyroot)
    # print?
    poles <- ts(poles, start = start(smooth), deltat = deltat(smooth) * step)
    poles
  })
}

print.ltv <- function(x, ...) {
  rx <- residuals(x)
  cat(paste("\nLinear transfer function model with \"", class(rx)[1],
    "\" data:\n",
    sep = ""
  ))
  cat(paste("Start = ", index2char(index(rx)[1], x$frequency),
    ", End = ", index2char(index(rx)[length(rx)], x$frequency),
    "\n",
    sep = ""
  ))
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  cat("Order: (na=", x$na, ", nb=", x$nb, ", nc=", x$nc,
    if (x$hasOffset) ", +offset" else ", no offset", ")\n",
    sep = ""
  )
  cat("LTV Parameters:\n")
  theta.etc <- rbind(
    mean = coef(x), std.dev. = x$sd.coef,
    deparse.level = 0
  )
  print.default(theta.etc, ..., print.gap = 2)
  invisible(x)
}

summary.ltv <- function(object, ...) {
  z <- object
  ans <- z[c("call", "P")]

  # diagnostics
  ans$msv <- mean(z$innov^2)
  ans$mse <- mean(z$resid^2)
  ans$v_e_ratio <- ans$msv / ans$mse

  vm <- z$y - z$sim.y # simulation-mode o/p. error
  ydiff <- diff(z$y) # for comparison with sim.-mode error

  # Mean-square simulation-mode output error
  ans$msvm <- mean(vm^2)
  ans$msydiff <- mean(ydiff^2)
  # m.s. simulation-mode error/m.s. change in flow
  ans$vm_ydiff_ratio <- ans$msvm / ans$msydiff

  # plot o/p. equation offset and d statistic
  # ts.plot(z$smooth$d, z$dstat)
  # title('change statistic and smoothed parameters,offset')

  class(ans) <- "summary.ltv"
  ans
}

print.summary.ltv <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )
  cat(
    "M.s. smoothed innovation:", formatC(x$msv, digits = digits),
    "\nM.s. smoothed residual:", formatC(x$mse, digits = digits),
    "\nRatio (M.s.v / M.s.e):", formatC(x$v_e_ratio, digits = digits),
    "\n\n"
  )
  cat(
    "Mean-square simulation-mode output error:", formatC(x$msvm, digits = digits),
    "\nMean-square change in output (ydiff):", formatC(x$msydiff, digits = digits),
    "\nRatio (M.s.vm / M.s.ydiff):", formatC(x$vm_ydiff_ratio, digits = digits),
    "\n\n"
  )
  cat("Final estimated error covariance matrix:\n")
  printCoefmat(x$P, digits = digits, ...)
  invisible(x)
}

vcov.ltv <- function(object, ...) {
  object$P
}

tsdiag.ltv <- function(object, gof.lag = 10, ...) {
  stats:::tsdiag.Arima(object, gof.lag, ...)
}
