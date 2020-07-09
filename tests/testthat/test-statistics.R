library(testthat)
library(hydromad)

context("Statistics and summary() methods")

data(SalmonBrook)
obsdat <- window(SalmonBrook, start = "1990-01-01", end = "1992-01-01")

## a fixed model to evaluate
mod <- hydromad(obsdat,
  sma = "cwi", tw = 30, f = 0.5,
  routing = "expuh", tau_s = 30, tau_q = 2, v_s = 0.3
)

test_that("basic summary() works", {
  s <- summary(mod, with.hydrostats = TRUE)
  expect_is(s, "summary.hydromad")
  expect_true(is.finite(s$rel.bias))
  #expect_equal(as.character(print(s)), "Time steps:") # is this test important?
  s.breaks <- summary(mod, breaks = "12 months")
  expect_is(s.breaks, "zoo")
})

test_that("all statistics can be evaluated", {
  hydromad_stats <- hydromad.stats()
  # Ignore r.sq.vartd because it requires event to be specified
  hydromad_stats$r.sq.vartd <- NULL
  ss <- objFunVal(mod, hydromad_stats)
  expect_is(ss, "list")
  ok <- sapply(ss, is.finite)
  if (!all(ok)) {
    warning(
      "non finite stats: ",
      toString(names(ss)[!ok])
    )
  }
  expect_true(all(is.finite(unlist(ss))))
})

test_that("custom objective functions work", {
  spec <- update(mod, v_s = c(0, 1))
  set.seed(0)
  fit1 <- fitByOptim1(spec, function(Q, X, ...) nseStat(Q, X))
  expect_is(fit1, "hydromad")
  set.seed(0)
  fit2 <- fitByOptim1(spec, ~ hmadstat("r.squared")(Q, X) - 1)
  expect_is(fit2, "hydromad")
  expect_equal(objFunVal(fit1), objFunVal(fit2))
  set.seed(0)
  fit3 <- fitByOptim1(spec, function(Q, X, ...) {
    hmadstat("r.sq.log")(Q, X) - 0.5 * hmadstat("rel.bias")(Q, X)
  })
  expect_is(fit3, "hydromad")
  expect_true(coef(fit1)[["v_s"]] != coef(fit3)[["v_s"]])
})

test_that("formula works within functions", {
  library(hydromad)
  data(HydroTestData)
  modx <- hydromad(HydroTestData, sma = "cmd", routing = "expuh", d = 200, f = 0.5, e = 0.1, tau_s = 10)

  # Formula gives correct answer outside functions
  expect_equal(objFunVal(modx, hmadstat("r.squared")), objFunVal(modx, ~ hmadstat("r.squared")(Q, X)))

  # Formula acknowledges collisions but overrides them with local values
  X <- "X should not cause an error"
  Q <- "Q should not cause an error"
  DATA <- "DATA should not cause an error"
  model <- "model should not cause an error"
  expect_warning(objFunVal(modx, ~ hmadstat("r.squared")(Q, X)))
  expect_equal(objFunVal(modx, hmadstat("r.squared")), objFunVal(modx, ~ hmadstat("r.squared")(Q, X)))

  formula.in.function <- function(w = 1) objFunVal(modx, ~ w * hmadstat("r.squared")(Q, X, ...))

  expect_error(formula.in.function(), NA)
  expect_equal(formula.in.function(), objFunVal(modx, hmadstat("r.squared")))

  # Formula is less than 10% slower
  time.fun <- system.time(replicate(1e3, objFunVal(modx, hmadstat("r.squared"))))
  time.formula <- system.time(replicate(1e3, objFunVal(modx, ~ hmadstat("r.squared")(Q, X, ...))))
  expect_lt((time.formula[3] - time.fun[3]) / time.fun[3], 0.1)
})
