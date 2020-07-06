library(testthat)
library(hydromad)

context("Soil Moisture Accounting models")

set.seed(0)
warmup <- 100
P <- ts(pmax(rnorm(200), 0))
E <- ts(20 + 10 * cos((1:200) / 20))
Q <- P * E * runif(P)
DATA <- cbind(P = P, E = E)
## some SMAs require Q also in simulation
DATAQ <- cbind(P = P, Q = Q)

test_that("cmd simulation is the same in R and C", {
  set.seed(0)
  mod0 <- hydromad(DATA, sma = "cmd")
  for (mod in simulate(mod0, 5)) {
    Csim <- predict(mod)
    expect_that(all(Csim >= 0), is_true())
    hydromad.options(pure.R.code = TRUE)
    pureRsim <- predict(mod)
    hydromad.options(pure.R.code = FALSE)
    expect_that(Csim, equals(pureRsim))
  }
})

test_that("cmd power-form simulation is the same in R and C", {
  set.seed(0)
  mod0 <- hydromad(DATA, sma = "cmd", shape = 2)
  for (mod in simulate(mod0, 5)) {
    Csim <- predict(mod)
    expect_that(all(Csim >= 0), is_true())
    hydromad.options(pure.R.code = TRUE)
    pureRsim <- predict(mod)
    hydromad.options(pure.R.code = FALSE)
    expect_that(Csim, equals(pureRsim))
  }
})

test_that("cwi simulation is the same in R and C", {
  set.seed(0)
  mod0 <- hydromad(DATA, sma = "cwi")
  for (mod in simulate(mod0, 5)) {
    Csim <- predict(mod)
    expect_that(all(Csim >= 0), is_true())
    hydromad.options(pure.R.code = TRUE)
    pureRsim <- predict(mod)
    hydromad.options(pure.R.code = FALSE)
    expect_that(Csim, equals(pureRsim))
  }
})

test_that("bucket simulation is the same in R and C", {
  set.seed(0)
  mod0 <- hydromad(DATA, sma = "bucket")
  for (mod in simulate(mod0, 5)) {
    Csim <- predict(mod)
    expect_that(all(Csim >= 0), is_true())
    hydromad.options(pure.R.code = TRUE)
    pureRsim <- predict(mod)
    hydromad.options(pure.R.code = FALSE)
    expect_that(Csim, equals(pureRsim))
  }
})

test_that("gr4j simulation is the same in R and C", {
  set.seed(0)
  mod0 <- hydromad(DATA,
    sma = "gr4j", routing = "gr4jrouting",
    etmult = 0.1, S_0 = 0.5, R_0 = 0
  )
  for (mod in simulate(mod0, 5)) {
    Csim <- na.trim(predict(mod))
    expect_that(all(Csim >= 0), is_true())
    hydromad.options(pure.R.code = TRUE)
    pureRsim <- na.trim(predict(mod))
    hydromad.options(pure.R.code = FALSE)
    expect_that(Csim, equals(pureRsim))
  }
})

test_that("AWBM simulation is the same in R and C", {
  set.seed(0)
  mod0 <- hydromad(DATA, sma = "awbm")
  for (mod in simulate(mod0, 5)) {
    Csim <- predict(mod)
    expect_that(all(Csim >= 0), is_true())
    hydromad.options(pure.R.code = TRUE)
    pureRsim <- predict(mod)
    hydromad.options(pure.R.code = FALSE)
    expect_that(Csim, equals(pureRsim))
  }
})

test_that("snow simulation is the same in R and C", {
  set.seed(0)
  DATA[, "E"] <- scale(DATA[, "E"]) * 2
  mod0 <- hydromad(DATA, sma = "snow")
  for (mod in simulate(mod0, 5)) {
    Csim <- predict(mod)
    expect_that(all(Csim >= 0), is_true())
    hydromad.options(pure.R.code = TRUE)
    pureRsim <- predict(mod)
    hydromad.options(pure.R.code = FALSE)
    expect_that(Csim, equals(pureRsim))
  }
})

test_that("simhyd simulation is the same in R and C", {
  set.seed(0)
  mod0 <- hydromad(DATA, sma = "simhyd")
  for (mod in simulate(mod0, 5)) {
    Csim <- predict(mod)
    expect_that(all(Csim >= 0), is_true())
    hydromad.options(pure.R.code = TRUE)
    pureRsim <- predict(mod)
    hydromad.options(pure.R.code = FALSE)
    expect_that(Csim, equals(pureRsim))
  }
})

test_that("sacramento simulation is the same in R and C", {
  set.seed(0)
  mod0 <- hydromad(DATA, sma = "sacramento")
  for (mod in simulate(mod0, 5)) {
    Csim <- predict(mod)
    expect_that(all(Csim >= 0), is_true())
    hydromad.options(pure.R.code = TRUE)
    pureRsim <- predict(mod)
    hydromad.options(pure.R.code = FALSE)
    expect_that(Csim, equals(pureRsim))
  }
})

test_that("scalar simulation runs", {
  set.seed(0)
  mod0 <- hydromad(DATA, sma = "scalar")
  for (mod in simulate(mod0, 5)) {
    expect_that(all(predict(mod) >= 0), is_true())
  }
})

test_that("runoffratio simulation runs", {
  set.seed(0)
  mod0 <- hydromad(DATAQ, sma = "runoffratio")
  for (mod in simulate(mod0, 5)) {
    expect_that(all(predict(mod) >= 0), is_true())
  }
})

test_that("dbm simulation runs", {
  set.seed(0)
  mod0 <- hydromad(DATAQ, sma = "dbm")
  for (mod in simulate(mod0, 5)) {
    expect_that(all(na.trim(predict(mod)) >= 0), is_true())
  }
})
