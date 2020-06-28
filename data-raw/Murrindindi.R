Murrindindi <- local({
  dat <- read.table("Murrindindi.dat", header = TRUE, as.is = TRUE)
  dat$Date <- as.Date(dat$Date)

  ## convert from ML/day to mm/day (i.e. divide by catchment area km^2)
  dat$Q <- dat$Q / 104.9
  #    dat$Q <- convertFlow(dat$Q, from = "ML", area.km2 = 104.9)

  ## zoo objects
  library(zoo)
  tsPQE <- zoo(dat[, -1], order.by = dat$Date, frequency = 1)
  tsPQE
})
usethis::use_data(Murrindindi, overwrite = TRUE)
