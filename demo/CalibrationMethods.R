library(hydromad)

data(SalmonBrook)

obsdat <- window(SalmonBrook, start = "1990-01-01", end = "1992-01-01")

P <- obsdat$P
simU <- hydromad.sim(obsdat,
  routing = NULL,
  sma = "cwi", tw = 30, f = 0.5, scale = 1 / 1000
)
simQ <- hydromad.sim(simU,
  sma = NULL,
  routing = "expuh", tau_s = 30, tau_q = 2, v_s = 0.3
)

simQ <- hydromad.sim(obsdat,
  sma = "cwi", tw = 30, f = 0.5, scale = 1 / 1000,
  routing = "expuh", tau_s = 30, tau_q = 2, v_s = 0.3
)


modDat <- merge(obsdat[, c("P", "E")], Q = byDays(simQ))
spec0 <- hydromad(modDat, sma = "cwi", routing = "uh")
rspecs <-
  list(
    ls21 = update(spec0, rfit = list("ls", order = c(2, 1))),
    sriv21 = update(spec0, rfit = list("sriv", order = c(2, 1))),
    inverse = update(spec0, rfit = list("inverse", order = c(2, 1)))
  )
#         sriv11 = update(spec0, rfit = list("sriv", order = c(1, 1))),
#         sriv20 = update(spec0, rfit = list("sriv", order = c(2, 0))),
#         sriv32 = update(spec0, rfit = list("sriv", order = c(3, 2))))
# cmdspecs <- lapply(cwispecs, update, sma = "cmd")

## first just test routing fitting with exact U given
fits <- as.runlist(lapply(rspecs, update, tw = 30, f = 0.5, scale = 1 / 1000))
## now test fitting of all parameters (SMA + expuh routing)
fullspec <- update(spec0, routing = "expuh", tau_q = c(0, 3), tau_s = c(3, 100), v_s = c(0, 1))
fits$Optim <- fitByOptim(fullspec)
fits$OptimBFGS <- fitByOptim(fullspec, method = "BFGS")
fits$Sampling <- fitBySampling(fullspec)
fits$SCE <- fitBySCE(fullspec)
fits$DE <- fitByDE(fullspec)

## OK, all fit correctly

## now test robustness of methods to:
## 1. error in input P (carries through to U via [correct] SMA)
## 2. error in SMA parameters or model structure
## 3. error in transfer function model order
## 4. error in streamflow Q


## objective function surface over parameters
sampLS21 <- simulate(rspecs$ls21, 200, sampletype = "random", FUN = objFunVal, bind = TRUE)
levelplot(result ~ tw * f, sampLS21, cex = 2, panel = panel.levelplot.points)

## NOTE: sampling is random, so results can vary!
set.seed(0)

# Q <- simulate(specs$cwiExpUH, nsim=1, sampletype = "random", FUN = predict)

fits <- list()
for (i in 1:3) {
  fits[[i]] <-
    with(
      specs,
      runlist(
        latinLS = fitBySampling(cwiLS21, samples = 100),
        latinSRIV = fitBySampling(cwiSRIV21, samples = 100),
        bfgsLS = fitByOptim(cwiLS21, samples = 64),
        bfgsSRIV = fitByOptim(cwiSRIV21, samples = 64),
        simplexLS = fitByOptim(specs$cwiLS21, samples = 64, method = "Nelder-Mead"),
        simplexSRIV = fitByOptim(specs$cwiSRIV21, samples = 64, method = "Nelder-Mead"),
        SCE = fitBySCE(cwiExpUH),
        DE = fitByDE(cwiExpUH)
      )
    )
}

lapply(fits, summary)

lapply(fits, function(x) {
  cbind(
    summary(x),
    t(sapply(x, function(x) {
      c(x["funevals"],
        time = x$timing[[1]]
      )
    }))
  )
})

lapply(fits, function(fit) {
  sapply(fit, function(x) c(x["funevals"], x$timing["elapsed"]))
})

## TODO: CMD model
## TODO: different simulation parameters
##       - third order routing
##       -
## TODO: noise model
## TODO: tf.invese.fit


modSpec <- hydromad(modDat,
  sma = "cwi",
  routing = "expuh", rfit = list("sriv", order = c(2, 1))
)

hydromad.options(trace = TRUE)

bfgsFit <- fitByOptim(modSpec, samples = 64, method = "BFGS")
summary(bfgsFit)
