library(testthat)
library(hydromad)
source("helper-ts-equals.R")
source("helper-expuh-orders.R")

context("Calibration methods on exact simulated data")

data(SalmonBrook)
obsdat <- window(SalmonBrook, start = "1990-01-01", end = "1992-01-01")

## joint simulation of SMA and routing -- the "true" model for testing
simQ <- fitted(hydromad(obsdat,
  sma = "cwi", tw = 30, f = 0.5, scale = 1 / 1000,
  routing = "expuh", tau_s = 30, tau_q = 2, v_s = 0.3
),
all = TRUE
)

modDat <- merge(obsdat[, c("P", "E")], Q = simQ)
spec <- hydromad(modDat, sma = "cwi", routing = "armax")
## give actual SMA parameter values, so we are fitting Q from exact U
xspec <- update(spec, tw = 30, f = 0.5, scale = 1 / 1000)


test_that("rfit methods work on (2,1) model with exact inputs", {
  fits <-
    runlist(
      ls = update(xspec, rfit = list("ls", order = c(2, 1))),
      sriv = update(xspec, rfit = list("sriv", order = c(2, 1))),
      inv = update(xspec, rfit = list("inverse", order = c(2, 1)))
    )
  expect_that(summary(fits$ls)$r.squared > 0.9999, is_true())
  expect_that(summary(fits$sriv)$r.squared > 0.9999, is_true())
  expect_that(summary(fits$inv)$r.squared > 0.9999, is_true())
})

U <- fitted(xspec, U = TRUE, all = TRUE)

Q_n0m0 <- expuh.sim(U, pars = n0m0)
Q_n1m0 <- expuh.sim(U, pars = n1m0)
Q_n1m1 <- expuh.sim(U, pars = n1m1)
Q_n2m0 <- expuh.sim(U, pars = n2m0)
Q_n2m1 <- expuh.sim(U, pars = n2m1)
Q_n2m2 <- expuh.sim(U, pars = n2m2)
Q_n3s0 <- expuh.sim(U, pars = n3s0)
Q_n3s1 <- expuh.sim(U, pars = n3s1)
Q_n3s2 <- expuh.sim(U, pars = n3s2)
Q_n3s3 <- expuh.sim(U, pars = n3s3)

test_that("Least squares (armax) fitting works for all orders with exact inputs", {
  ## note hydromad.getOption("warmup") == 100
  expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n0m0), order = c(0, 0))), ts_equals(Q_n0m0, tol = 1e-5))
  expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n1m0), order = c(1, 0))), ts_equals(Q_n1m0, tol = 1e-5))
  expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n1m1), order = c(1, 1))), ts_equals(Q_n1m1, tol = 1e-5))
  expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n2m0), order = c(2, 0))), ts_equals(Q_n2m0, tol = 1e-5))
  expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n2m1), order = c(2, 1))), ts_equals(Q_n2m1, tol = 1e-5))
  expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n2m2), order = c(2, 2))), ts_equals(Q_n2m2, tol = 1e-5))
  expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n3s0), order = c(3, 2))), ts_equals(Q_n3s0, tol = 1e-5))
  expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n3s1), order = c(3, 2))), ts_equals(Q_n3s1, tol = 1e-5))
  expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n3s2), order = c(3, 1))), ts_equals(Q_n3s2, tol = 1e-5))
  expect_that(fitted(armax.ls.fit(cbind(U = U, Q = Q_n3s3), order = c(3, 0), delay = 0)), ts_equals(Q_n3s3, tol = 1e-5))
})

test_that("SRIV (armax) fitting works for all orders with exact inputs", {
  ## note hydromad.getOption("warmup") == 100

  expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n0m0), order = c(0, 0))), ts_equals(Q_n0m0, tol = 1e-5))
  expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n1m0), order = c(1, 0))), ts_equals(Q_n1m0, tol = 1e-5))
  expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n1m1), order = c(1, 1))), ts_equals(Q_n1m1, tol = 1e-5))
  expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n2m0), order = c(2, 0))), ts_equals(Q_n2m0, tol = 1e-5))
  expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n2m1), order = c(2, 1))), ts_equals(Q_n2m1, tol = 1e-5))
  expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n2m2), order = c(2, 2))), ts_equals(Q_n2m2, tol = 1e-5))
  expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n3s0), order = c(3, 2))), ts_equals(Q_n3s0, tol = 1e-5))
  expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n3s1), order = c(3, 2))), ts_equals(Q_n3s1, tol = 1e-5))
  expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n3s2), order = c(3, 1))), ts_equals(Q_n3s2, tol = 1e-5))
  expect_that(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n3s3), order = c(3, 0), delay = 0)), ts_equals(Q_n3s3, tol = 1e-5))

  ## tmp <- ts.intersect(fitted(armax.sriv.fit(cbind(U = U, Q = Q_n3s3), order = c(3,0), delay = 0)), Q_n3s3)
  ## summary(tmp[,1] - tmp[,2])
  ## xyplot(ts(tmp), superpose = TRUE)
})

test_that("Inverse (armax) fitting works for all orders with exact inputs", {
  ## TODO: test use.Qm = TRUE; other options?
  hydromad.options(inverse.rel.tolerance = 1e-4)
  expect_that(fitted(armax.inverse.fit(cbind(P = obsdat$P, Q = Q_n0m0), order = c(0, 0), use.Qm = FALSE)), ts_equals(Q_n0m0, tol = 1e-5))
  expect_that(fitted(armax.inverse.fit(cbind(P = obsdat$P, Q = Q_n1m0), order = c(1, 0), use.Qm = FALSE)), ts_equals(Q_n1m0, tol = 1e-2))
  expect_that(fitted(armax.inverse.fit(cbind(P = obsdat$P, Q = Q_n1m1), order = c(1, 1), use.Qm = FALSE)), ts_equals(Q_n1m1, tol = 1e-2))
  expect_that(fitted(armax.inverse.fit(cbind(P = obsdat$P, Q = Q_n2m0), order = c(2, 0), use.Qm = FALSE)), ts_equals(Q_n2m0, tol = 1e-2))
  expect_that(fitted(armax.inverse.fit(cbind(P = obsdat$P, Q = Q_n2m1), order = c(2, 1), use.Qm = FALSE)), ts_equals(Q_n2m1, tol = 1e-2))
  expect_that(fitted(armax.inverse.fit(cbind(P = obsdat$P, Q = Q_n2m2), order = c(2, 2), use.Qm = FALSE)), ts_equals(Q_n2m2, tol = 1e-2))
  expect_that(fitted(armax.inverse.fit(cbind(P = obsdat$P, Q = Q_n3s0), order = c(3, 2), use.Qm = FALSE)), ts_equals(Q_n3s0, tol = 1e-2))
  expect_that(fitted(armax.inverse.fit(cbind(P = obsdat$P, Q = Q_n3s1), order = c(3, 2), use.Qm = FALSE)), ts_equals(Q_n3s1, tol = 1e-2))
  expect_that(fitted(armax.inverse.fit(cbind(P = obsdat$P, Q = Q_n3s2), order = c(3, 1), use.Qm = FALSE)), ts_equals(Q_n3s2, tol = 1e-2))
  expect_that(fitted(armax.inverse.fit(cbind(P = obsdat$P, Q = Q_n3s3), order = c(3, 0), delay = 0, use.Qm = FALSE)), ts_equals(Q_n3s3, tol = 1e-2))
  ## tmp <- cbind(fitted(armax.inverse.fit(cbind(P = obsdat$P, Q = Q_n3s3), order = c(3,0), delay = 0, use.Qm = FALSE)), Q_n3s3)
  ## summary(tmp[,1] - tmp[,2])
  ## xyplot(ts(tmp), superpose = TRUE)
})

test_that("SMA joint fitting methods work with exact inputs", {
  jspec <- update(spec,
    routing = "expuh",
    tau_q = c(0, 3), tau_s = c(3, 100), v_s = c(0, 1)
  )
  ## TODO: suppress optim() output
  set.seed(0)
  expect_that(summary(fitByOptim(jspec, control = list(maxit = 150)))$r.squared > 0.98, is_true())
  set.seed(0)
  expect_that(summary(fitBySCE(jspec, control = list(maxeval = 150)))$r.squared > 0.98, is_true())
  set.seed(0)
  expect_that(summary(fitByDE(jspec, control = list(itermax = 4)))$r.squared > 0.98, is_true())
  set.seed(0)
  expect_that(summary(fitByDream(jspec, control = list(ndraw = 150)))$r.squared > 0.98, is_true())
  set.seed(0)
  # TODO: why is CMAES performing less well?
  expect_that(summary(fitByCMAES(jspec, control = list(maxit = 20)))$r.squared > 0.96, is_true())
  set.seed(0)
  expect_that(summary(fitByDDS(jspec, control = modifyList(hydromad.getOption("dds.control"), list(max_number_function_calls = 150))))$r.squared > 0.98, is_true())
  set.seed(0)
  expect_that(summary(fitByNsga2(jspec, control = list(generations = 2)))$r.squared > 0.98, is_true())
})
