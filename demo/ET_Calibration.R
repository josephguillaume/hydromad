################################################################################
## Prepare example data
data(Cotter)
x <- Cotter[1:1000]
## Define exact model to be used
fitt <- hydromad(
  DATA = x, sma = "cmd", f = 0.6, e = 0.166, d = 220, return_state = T,
  routing = "expuh", tau_q = 2, tau_s = 25, v_s = 0.1
)
## Calculate observed actual ET
x$aET <- fitt$U$ET
## Calculate exact Q resulting from calculated ET
x$Q <- fitted(fitt, all = TRUE)

################################################################################
## Method 1.
## First fit SMA using actual ET
mod.sma <- hydromad(DATA = x, sma = "cmd", return_state = T)
fit.sma <- fitByOptim(mod.sma, objective = ~ hmadstat("r.squared")(DATA$aET, U$ET))
fit.sma

## Then fit routing using Q
mod <- update(fit.sma, routing = "expuh", tau_q = c(0, 10), tau_s = c(10, 50), v_s = c(0, 1))
fit <- fitByOptim(mod, objective = hmadstat("r.squared"))
fit

################################################################################
## Method 2.
## Fit both SMA and routing using Q
##  with e fixed to improve identifiability
modo <- hydromad(
  DATA = x, sma = "cmd", e = 0.166,
  routing = "expuh", tau_q = c(0, 10), tau_s = c(10, 50), v_s = c(0, 1)
)
hydromad.options(objective = hmadstat("r.squared"))
fito <- fitByOptim(modo)
fito

################################################################################
## Method 3.
## Fit SMA using actual ET, and fit routing using
##   SRIV with inverse simulation from the observed streamflow data
modr <- hydromad(
  DATA = x, sma = "cmd", return_state = T,
  routing = "expuh", rfit = list("inverse", order = c(2, 1))
)
fitr <- fitByOptim(modr, objective = ~ hmadstat("r.squared")(DATA$aET, U$ET))
fitr

## results
summary(runlist(fit, fito, fitr))
