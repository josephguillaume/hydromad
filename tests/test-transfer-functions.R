library(testthat)
library(hydromad)
source("helper-expuh-orders.R")

context("Transfer Functions")

test_that("cwi+expuh(2,1) simulation looks reasonable", {
  data(SalmonBrook)
  obsdat <- window(SalmonBrook, start = "1990-01-01", end = "1992-01-01")
  ## joint simulation of SMA and routing
  simQ <-
    fitted(hydromad(obsdat,
      sma = "cwi", tw = 30, f = 0.5, scale = 1 / 1000,
      routing = "expuh", tau_s = 30, tau_q = 2, v_s = 0.3
    ),
    all = TRUE
    )
  expect_that(NROW(simQ) == NROW(obsdat), is_true())
  expect_that(all(is.finite(simQ)), is_true())
})

test_that("TF parameter conversions are consistent", {
  toAbAndBack <- function(expars) {
    tfParsConvert(tfParsConvert(expars, "a,b"), "tau,v")[names(expars)]
  }
  expect_that(toAbAndBack(n0m0), equals(n0m0))
  expect_that(toAbAndBack(n1m0), equals(n1m0))
  expect_that(toAbAndBack(n1m1), equals(n1m1))
  expect_that(toAbAndBack(n2m0), equals(n2m0))
  expect_that(toAbAndBack(n2m1), equals(n2m1))
  expect_that(toAbAndBack(n2m2), equals(n2m2))
  expect_that(toAbAndBack(n3s0), equals(n3s0))
  # expect_that(toAbAndBack(n3s1), equals(n3s1)) ## can't deduce series=1|0
  expect_that(toAbAndBack(n3s2), equals(n3s2))
  expect_that(toAbAndBack(n3s3), equals(n3s3))

  set.seed(0)
  warmup <- 100
  U <- ts(pmax(rnorm(200), 0))
  armaxSim <- function(expars, ..., init = 0) {
    window(armax.sim(U, pars = tfParsConvert(expars, "a,b"), ..., init = init),
      start = warmup, end = 198
    )
  }
  expuhSim <- function(expars, ...) {
    window(expuh.sim(U, pars = expars, ...),
      start = warmup, end = 198
    )
  }

  ## there are small(ish) differences due to initialisation
  tol <- 0.05
  expect_that(armaxSim(n0m0), equals(expuhSim(n0m0), tol = tol))
  expect_that(armaxSim(n1m0), equals(expuhSim(n1m0), tol = tol))
  expect_that(armaxSim(n1m1), equals(expuhSim(n1m1), tol = tol))
  expect_that(armaxSim(n2m0), equals(expuhSim(n2m0), tol = tol))
  expect_that(armaxSim(n2m1), equals(expuhSim(n2m1), tol = tol))
  expect_that(armaxSim(n2m2), equals(expuhSim(n2m2), tol = tol))
  expect_that(armaxSim(n3s0), equals(expuhSim(n3s0), tol = tol))
  expect_that(armaxSim(n3s1), equals(expuhSim(n3s1), tol = tol))
  expect_that(armaxSim(n3s2), equals(expuhSim(n3s2), tol = tol))
  expect_that(armaxSim(n3s3), equals(expuhSim(n3s3), tol = tol))
  # xyplot(cbind(armaxSim(n2m2), expuhSim(n2m2)), superpose = TRUE)
})


test_that("TF simulation and inverse simulation are consistent", {
  set.seed(0)
  warmup <- 100
  U <- ts(pmax(rnorm(200), 0))
  P <- U + abs(rnorm(U)) / 3
  Uw <- window(U, start = warmup, end = end(U)[1] - 2)
  Pw <- window(P, start = warmup, end = end(U)[1] - 2)
  expuhSim <- function(expars, ...) {
    window(expuh.sim(U, pars = expars, ...),
      start = warmup, end = end(U)[1] - 2
    )
  }
  for (pure.R.code in c(TRUE, FALSE)) {
    message(paste("testing with pure.R.code =", pure.R.code))
    hydromad.options(pure.R.code = pure.R.code)
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n0m0), P = Pw), pars = tfParsConvert(n0m0), use.Qm = FALSE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n1m0), P = Pw), pars = tfParsConvert(n1m0), use.Qm = FALSE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n1m1), P = Pw), pars = tfParsConvert(n1m1), use.Qm = FALSE)[-(1:3)], equals(Uw[-(1:3)]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n2m0), P = Pw), pars = tfParsConvert(n2m0), use.Qm = FALSE)[-(1:3)], equals(Uw[-(1:3)]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n2m1), P = Pw), pars = tfParsConvert(n2m1), use.Qm = FALSE)[-(1:3)], equals(Uw[-(1:3)]))
  }
  ## third order models, and all simulations with use.Qm = TRUE, seem to need a longer warmup
  set.seed(0)
  warmup <- 200
  U <- ts(pmax(rnorm(300), 0))
  P <- U + abs(rnorm(U)) / 3
  Uw <- window(U, start = warmup, end = end(U)[1] - 2)
  Pw <- window(P, start = warmup, end = end(U)[1] - 2)
  for (pure.R.code in c(TRUE, FALSE)) {
    message(paste("testing with pure.R.code =", pure.R.code))
    hydromad.options(pure.R.code = pure.R.code)
    ## higher order models with use.Qm = FALSE (completing the previous section)
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n2m2), P = Pw), pars = tfParsConvert(n2m2), use.Qm = FALSE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n3s0), P = Pw), pars = tfParsConvert(n3s0), use.Qm = FALSE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n3s1), P = Pw), pars = tfParsConvert(n3s1), use.Qm = FALSE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n3s2), P = Pw), pars = tfParsConvert(n3s2), use.Qm = FALSE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n3s3), P = Pw), pars = tfParsConvert(n3s3), use.Qm = FALSE)[-1], equals(Uw[-1]))
    ## all models with use.Qm = TRUE
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n0m0), P = Pw), pars = tfParsConvert(n0m0), use.Qm = TRUE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n1m0), P = Pw), pars = tfParsConvert(n1m0), use.Qm = TRUE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n1m1), P = Pw), pars = tfParsConvert(n1m1), use.Qm = TRUE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n2m0), P = Pw), pars = tfParsConvert(n2m0), use.Qm = TRUE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n2m1), P = Pw), pars = tfParsConvert(n2m1), use.Qm = TRUE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n2m2), P = Pw), pars = tfParsConvert(n2m2), use.Qm = TRUE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n3s0), P = Pw), pars = tfParsConvert(n3s0), use.Qm = TRUE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n3s1), P = Pw), pars = tfParsConvert(n3s1), use.Qm = TRUE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n3s2), P = Pw), pars = tfParsConvert(n3s2), use.Qm = TRUE)[-1], equals(Uw[-1]))
    expect_that(armax.inverse.sim(cbind(Q = expuhSim(n3s3), P = Pw), pars = tfParsConvert(n3s3), use.Qm = TRUE)[-1], equals(Uw[-1]))
  }
  # armax.inverse.sim(cbind(Q = expuhSim(n1m1), P = Pw), pars = tfParsConvert(n1m1)) - Uw
  # armax.inverse.sim(cbind(Q = expuhSim(n3s0), P = Pw), pars = tfParsConvert(n3s0)) - Uw
  # xyplot(cbind(armax.inverse.sim(cbind(Q = expuhSim(n1m0), P = Pw), pars = tfParsConvert(n1m0)), Uw), superpose = TRUE)
  # xyplot(cbind(armax.inverse.sim(cbind(Q = expuhSim(n1m1), P = Pw), pars = tfParsConvert(n1m1)), Uw), superpose = TRUE)
  # xyplot(cbind(armax.inverse.sim(cbind(Q = expuhSim(n3s0), P = Pw), pars = tfParsConvert(n3s0)), Uw), superpose = TRUE)
})
