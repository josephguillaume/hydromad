
library(hydromad)
data(Queanbeyan)
ts83 <- window(Queanbeyan, start = "1983-01-01", end = "1986-01-01")

cmdspec80s <-
  hydromad(ts83,
    sma = "cmd",
    f = c(0.1, 1), e = 0.166, d = 200,
    shape = 2^(-1:5), ## use discrete values of shape on log scale, for simulation
    routing = "expuh", delay = 1,
    tau_s = c(5, 300), tau_q = c(0, 5), v_s = c(0, 1)
  )

## Calibrate a model specification with 7 methods, to show resulting optimisation trace plot.
tryAllFitMethods <-
  function(spec, maxeval, samples = hydromad.getOption("fit.samples")) {
    oopt <- hydromad.options(
      fit.samples = samples, optim.control = list(maxit = maxeval - samples),
      nlminb.control = list(iter.max = maxeval, eval.max = (maxeval - samples) / 5)
    )
    on.exit(hydromad.options(oopt))
    ## try all fitting methods with about 'maxeval' function evaluations
    meths <- list()
    message("__PORT__")
    set.seed(0)
    meths$PORT <- fitByOptim(spec, method = "PORT")
    message("__BFGS__")
    set.seed(0)
    meths$BFGS <- fitByOptim(spec, method = "BFGS")
    message("__NM__")
    set.seed(0)
    meths$NM <- fitByOptim(spec, method = "Nelder-Mead")
    message("__SANN__")
    set.seed(0)
    meths$SANN <- fitByOptim(spec, method = "SANN")
    message("__SCE__")
    set.seed(0)
    meths$SCE <- fitBySCE(spec, control = list(maxeval = maxeval))
    message("__DE__")
    set.seed(0)
    meths$DE <- fitByDE(spec, control = list(itermax = maxeval / 50))
    message("__DREAM__")
    set.seed(0)
    meths$DREAM <- fitByDream(spec, control = list(ndraw = maxeval))
    as.runlist(meths)
  }
hydromad.options(objective = hmadstat("r.squared"))
hydromad.options(loglik = ~ -0.5 * sum((Q - X)^2, na.rm = TRUE))
## three-store routing version, for comparing optimisation algorithms
cmd3spec <- update(cmdspec80s, tau_3 = c(1, 150), v_3 = c(0, 0.5))
cmd3fits <- tryAllFitMethods(cmd3spec, maxeval = 5000, samples = 100)
## calculate optimisation traces
everyN <- function(z, N = 50) {
  aggregate(z, function(x) N * ceiling(x / N),
    FUN = function(x) if (all(is.na(x))) NA else tail(na.omit(x), 1)
  )
}
traces <- lapply(cmd3fits, optimtrace, objective = hmadstat("r.squared"))
traces <- everyN(na.approx(do.call("merge", traces), na.rm = FALSE))
save(traces, file = "hydromad_optimtraces.Rdata")


## Latin Hypercube sampling of parameter space and calculating fit statistics
set.seed(0)
sims <- simulate(cmdspec80s, 2000,
  sampletype = "latin",
  objective = hydromad.stats(c("r.squared", "r.sq.log"))
)
## remove fixed parameters
sims <- sims[, -c(2:3)]
## calculate 90% coverage sets
gluer2 <- defineFeasibleSet(sims[, 1:6],
  model = cmdspec80s,
  objseq = sims$r.squared,
  target.coverage = 0.9
)
gluer2log <- defineFeasibleSet(sims[, 1:6],
  model = cmdspec80s,
  objseq = sims$r.sq.log,
  target.coverage = 0.9
)
save(sims, gluer2, gluer2log, file = "hydromad_sims.Rdata")


## two different model structures
fits83.spec <-
  list(
    cmd = cmdspec80s,
    intensity = update(cmdspec80s, sma = "intensity", maxP = c(20, 1000), scale = c(0, 1))
  )
## fit each model using NSE
set.seed(0)
fits83.raw <-
  as.runlist(lapply(fits83.spec, fitByOptim, objective = hmadstat("r.squared")))
## fit each model using NSE_log
set.seed(0)
fits83.log <-
  as.runlist(lapply(fits83.spec, fitByOptim, objective = hmadstat("r.sq.log")))
save(fits83.raw, fits83.log, file = "hydromad_fits83.Rdata")
