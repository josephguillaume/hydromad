YeAl97 <- local({
  library(reshape)

  stat <- read.csv("YeAl97.csv", stringsAsFactors = FALSE)
  stat <- melt(stat, id.vars = c("Catchment", "calib.period", "sim.period", "perf.stat"), variable_name = "Model.str")
  stat <- cast(stat, ... ~ perf.stat)
  stat$E[is.na(stat$E)] <- -Inf
  stat <- as.data.frame(stat)
  stat
})
usethis::use_data(YeAl97, overwrite = TRUE)
