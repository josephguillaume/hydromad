library(testthat)
library(hydromad)


context("Calibration results with real data")


test_that("cwi+expuh optim gives reasonable result", {
  data(Cotter)
  x <- Cotter[1:1000]
  modx <- hydromad(x,
    sma = "cwi", routing = "expuh",
    tau_s = c(1, 100), v_s = c(0, 1)
  )
  ## now try to fit it
  set.seed(0)
  fitx <- fitByOptim(modx)
  s <- summary(fitx)
  expect_gt(s$r.squared, 0.65)
  expect_gt(s$r.sq.log, 0.8)
  expect_lt(s$rel.bias, 0.01)
})


test_that("cwi+armax optim/SRIV gives reasonable result", {
  data(Cotter)
  x <- Cotter[1:1000]
  modx <- hydromad(x,
    sma = "cwi", routing = "armax",
    rfit = list("sriv", order = c(2, 1))
  )
  ## now try to fit it
  set.seed(0)
  fitx <- fitByOptim(modx)
  s <- summary(fitx)
  expect_gt(s$r.squared, 0.75)
  expect_gt(s$r.sq.log, 0.85)
  expect_lt(abs(s$rel.bias), 0.07)
})




test_that("cmd+powuh SCE gives reasonable result", {
  data(Cotter)
  x <- Cotter[1:1000]
  modx <- hydromad(x, sma = "cmd", routing = "powuh")
  ## now try to fit it
  set.seed(0)
  fitx <- fitBySCE(modx, control = list(maxit = 3))
  s <- summary(fitx)
  expect_gt(s$r.squared, 0.5)
  expect_gt(s$r.sq.log, 0.6)
  expect_lt(abs(s$rel.bias), 0.06)
})
