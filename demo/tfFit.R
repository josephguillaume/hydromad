
## Simulation study of results from SRIV
## (Simple Refined Instrumental Variable) and
## Least Squares methods for fitting transfer functions.

## Reproduce results from Taylor et. al. (2007).
## Environmental time series analysis and forecasting
## with the Captain Toolbox.
## Environmental Modelling and Software 22, p. 811.

## A first-order autoregressive system, with two unit impulse inputs,
## with a serially correlated noise model (an AR(2) model).

library(hydromad)

U <- c(rep(0, 10), 1, rep(0, 50), 1, rep(0, 50)) ## Input signal
X <- filter(2 * U, 0.8, "recursive") ## Noise free output

ivCoef <- lsCoef <- invCoef <-
  matrix(NA, nrow = 1000, ncol = 2)

hydromad.options(prefilter = 0.5)

pb <- txtProgressBar(max = 1000, style = 3)
for (i in 1:1000) {
  E <- rnorm(length(U)) * 0.02 ## White noise
  Y <- X + filter(E, c(1.2, -0.8), "recursive") ## Stochastic output
  est1 <- tf.ls.fit(cbind(U = U, Q = Y),
    normalise = FALSE,
    order = c(n = 1, m = 0), delay = 0, warmup = 0
  ) ## LS estimation
  est2 <- tf.sriv.fit(cbind(U = U, Q = Y),
    normalise = FALSE,
    order = c(n = 1, m = 0), delay = 0, warmup = 0, epsilon = 1e-8
  ) ## SRIV estimation
  lsCoef[i, ] <- coef(est1, "tau,v")
  ivCoef[i, ] <- coef(est2, "tau,v")
  setTxtProgressBar(pb, i)
}
close(pb)

## True values are
c(10, -1 / log(0.8))

colMeans(ivCoef)
colMeans(lsCoef)

## show spread of fitted values around true values
tauhist <- histogram(~ ivCoef[, 1], nint = 25, sub = "using SRIV") +
  layer(panel.abline(v = 10, lwd = 3, lty = 3))
vhist <- histogram(~ ivCoef[, 2], nint = 25, sub = "using SRIV") +
  layer(panel.abline(v = -1 / log(0.8), lwd = 3, lty = 3))
print(c("Time Constant" = tauhist, "Steady State Gain" = vhist))

## compare spread of SRIV vs LS fits
taufits <- make.groups(SRIV = ivCoef[, 1], "Least Squares" = lsCoef[, 1])
vfits <- make.groups(SRIV = ivCoef[, 2], "Least Squares" = lsCoef[, 2])
taubwp <- bwplot(which ~ data, taufits)
vbwp <- bwplot(which ~ data, vfits)
print(c("Time Constant" = taubwp, "Steady State Gain" = vbwp))

## simulated response trace using mean parameter values,
## together with a sample realisation of the model.
tau_s <- mean(ivCoef[, 1])
v_s <- mean(ivCoef[, 2])
print(xyplot(Y) +
  layer(panel.lines(expuh.sim(U, tau_s = tau_s, v_s = v_s))))
