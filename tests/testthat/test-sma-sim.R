library(testthat)
library(hydromad)
library(patrick)

context("Soil Moisture Accounting models")

# Set up
set.seed(0)
warmup <- 100
P <- ts(pmax(rnorm(200), 0))
E <- ts(20 + 10 * cos((1:200) / 20))
Q <- P * E * runif(P)
Temp <- ts(10 + 15 * cos((1:200) / 20))
DATA <- cbind(P = P, E = E)
DATAE <- DATA
DATAE[, "E"] <- scale(DATAE[, "E"]) * 2
DATATE <- cbind(P = P, E = E * 0.5, T = Temp)
## some SMAs require Q also in simulation
DATAQ <- cbind(P = P, Q = Q)


with_parameters_test_that(
  "cmds are the same in R and C",
  {
    set.seed(0)
    n_iter <- 5
    for (mod in simulate(mod, n_iter)) {
      if (test_name == "sma = \"gr4j\"") {
        Csim <- na.trim(predict(mod))
        expect_true(all(Csim >= 0))
        hydromad.options(pure.R.code = TRUE)
        pureRsim <- na.trim(predict(mod))
        hydromad.options(pure.R.code = FALSE)
        expect_equal(Csim, pureRsim)
      }
      else {
        Csim <- predict(mod)
        expect_true(all(Csim >= 0))
        hydromad.options(pure.R.code = TRUE)
        pureRsim <- predict(mod)
        hydromad.options(pure.R.code = FALSE)
        expect_equal(Csim, pureRsim)
      }
    }
  },
  cases(
    `sma = "cmd"` = list(mod = hydromad(DATA, sma = "cmd")),
    `sma = "cmd", shape = 2` = list(mod = hydromad(DATA, sma = "cmd", shape = 2)),
    `sma = "cwi"` = list(mod = hydromad(DATA, sma = "cwi")),
    `sma = "bucket"` = list(mod = hydromad(DATA, sma = "bucket")),
    `sma = "awbm"` = list(mod = hydromad(DATA, sma = "awbm")),
    `sma = "gr4j"` = list(mod = hydromad(DATA,
      sma = "gr4j", routing = "gr4jrouting",
      etmult = 0.1, S_0 = 0.5, R_0 = 0
    )),
    `sma = "awbm"` = list(mod = hydromad(DATA, sma = "awbm")),
    `sma = "simhyd"` = list(mod = hydromad(DATA, sma = "simhyd")),
    `sma = "sacramento"` = list(mod = hydromad(DATA, sma = "sacramento")),
    `sma = "snow"` = list(mod = hydromad(DATAE, sma = "snow")),
    `sma = "hbv"` = list(mod = hydromad(DATATE,
      sma = "hbv", routing = "hbvrouting"
    ))
  )
)


with_parameters_test_that(
  "cmds run",
  {
    set.seed(0)
    n_iter <- 5
    for (mod in simulate(mod, n_iter)) {
      expect_true(all(na.trim(predict(mod)) >= 0))
    }
  },
  cases(
    `sma = "scalar"` = list(mod = hydromad(DATA, sma = "scalar")),
    `sma = "runoffratio"` = list(mod = hydromad(DATAQ, sma = "runoffratio")),
    `sma = "dbm"` = list(mod = hydromad(DATAQ, sma = "dbm"))
  )
)
