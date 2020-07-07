
library(hydromad)
data(Queanbeyan)

stop("this needs to be updated")

ts83 <- window(Queanbeyan, start = "1983-01-01", end = "1986-01-01")

cmdspec80s <-
  hydromad(ts83,
    sma = "cmd",
    f = c(0.1, 1),
    shape = 2^(-1:5),
    e = 0.166, d = 200,
    routing = "expuh",
    tau_s = c(5, 200), tau_q = c(0, 5), v_s = c(0, 1)
  )

sims <- simulate(cmdspec80s, 800,
  sampletype = "latin",
  objective = hydromad.stats(c(
    "r.squared", "r.sq.log"
  ))
)

## remove fixed parameters
sims <- sims[, -c(2:3)]

gluer2 <- glue.hydromad.default(cmdspec80s,
  objseq = sims$r.squared,
  psets = sims[, 1:5], target.coverage = 0.9
)
gluer2log <- glue.hydromad.default(cmdspec80s,
  objseq = sims$r.sq.log,
  psets = sims[, 1:5], target.coverage = 0.9
)

## show hydrographs with 90% coverage bounds from "glue" (not exactly GLUE method)

xyplot(gluer2, feasible.bounds = TRUE, with.P = TRUE)

xyplot(gluer2log, feasible.bounds = TRUE)

xyplot(gluer2log,
  feasible.bounds = TRUE, cut = 3,
  scales = list(y = list(log = TRUE), format = "%Y %b"),
  yscale.components = yscale.components.log10.3
)

## plot distributions of parameters in the 90% coverage sets

foo <- make.groups(
  r.squared = data.frame(gluer2$feasible.set),
  r.sq.log = data.frame(gluer2log$feasible.set),
  nonfeasible = sims[, 1:5]
)

qqmath(~ f + shape + tau_s + tau_q + v_s, foo,
  outer = TRUE,
  groups = which, auto.key = TRUE,
  scales = list(y = list(relation = "free"))
)

## density plot; this should ideally be weighted by the objective function values
# marginal.plot(foo[,1:5], groups = foo$which,
#              auto.key = list(title = "Statistic used to define\n90% coverage set", cex.title = 1,
#              corner = c(.9, 0.1)))

## pairwise scatterplots of parameter values, with a convex hull drawn around the "90% coverage set"

xyplot(shape + tau_s + tau_q + v_s ~ f, sims,
  outer = TRUE, scales = "free"
) +
  layer(subset <- sims$r.squared > 0.5, ii <- chull(x[subset], y[subset]), panel.polygon(x[subset][ii], y[subset][ii])) +
  layer(panel.levelplot.points(x, y, z = pmax(sims$r.squared, 0))) +
  layer(i <- which.max(sims$r.squared), panel.points(x[i], y[i], pch = 16, col = "black", cex = 2))

## all pairwise combinations

splom(foo[, 1:5],
  groups = foo$which,
  par.settings = simpleTheme(lwd = 2, pch = c(16:17)),
  auto.key = list(lines = TRUE),
  axis.text.cex = 0.7, jitter.x = TRUE, jitter.y = TRUE,
  panel = function(...) NULL,
  lower.panel = function(...) NULL
) +
  glayer(panel.xyplot(..., pch = ".", col.symbol = "grey"), groups = 3) +
  glayer(ii <- chull(x, y),
    panel.polygon(x[ii], y[ii], ..., col = NA, border = col.line),
    panel.xyplot(x[1:3], y[1:3], jitter.x = FALSE, jitter.y = FALSE, ...),
    groups = 1:2
  )
