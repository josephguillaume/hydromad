Wye <- local({
  dat <- read.table("Wye.dat",
    sep = "\t",
    header = TRUE, as.is = TRUE
  )
  dat$Time <- as.POSIXct(dat$Time, tz = "GMT")

  ## streamflow already in mm.
  ## for reference, catchment area is 10.6 km^2

  ## zoo objects
  library(zoo)
  tsPQ <- zoo(dat[, -1], order.by = dat$Time, frequency = 1)
  tsPQ
})
usethis::use_data(Wye, overwrite = TRUE)
