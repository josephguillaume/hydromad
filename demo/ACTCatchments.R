
library(hydromad)

## ACT CATCHMENTS
data(Cotter)
data(Queanbeyan)
data(Molonglo)
data(Orroral)
dsets <-
  c(
    "Cotter",
    "Queanbeyan",
    "Molonglo",
    "Orroral"
  )
names(dsets) <- dsets
dsets <- lapply(dsets, get)
## define calibration periods
calts <- lapply(dsets, window, start = "1970-01-01", end = "1980-01-01")
## calibrations
hydromad.options(order = c(n = 2, m = 1))
mods <- lapply(calts, function(DATA) {
  hydromad(DATA, rfit = "sriv")
})
mods <- as.runlist(mods)
summary(mods)
## simulations
smods <- mods
for (nm in names(mods)) {
  smods[[nm]] <- update(smods[[nm]], newdata = dsets[[nm]])
}
summary(smods, pars = FALSE)

## flow duration curves
xyplot.list(smods,
  FUN = qqmath, trans = log,
  f.value = ppoints(100), tails.n = 50,
  type = c("g", "b")
)
