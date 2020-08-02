library(testthat)
library(hydromad)

context("Smoke tests")

data(Cotter)
ts90s <- window(Cotter$Q, start = "1990-01-01", end = "1999-12-31")

## Test that tsdiag.hydromad returns correct class
test_that("can plot arima models with tsdiag", {
  expect_type(tsdiag(arima(Cotter$P)), "NULL")
  expect_type(tsdiag(arima(ts90s)), "NULL")
  expect_error(tsdiag(Cotter))
})

## Test that xyplot.hydromad works returns correct class
test_that("can plot data with xyplot plot", {
  expect_type(xyplot(Cotter), "trellis")
  expect_type(xyplot(ts90s), "trellis")
  expect_error(xyplot(arima(ts90s)))
})
