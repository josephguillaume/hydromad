Canning <- local({
  dat <- read.table("Canning.csv",
    sep = ",",
    col.names = c("Day", "Mon", "Year", "P", "Q", "E")
  )
  dat$Date <- with(dat, as.Date(ISOdate(Year, Mon, Day)))

  ## streamflow already in mm.
  ## for reference, catchment area is 517 km^2

  ## zoo objects
  library(zoo)
  tsPQE <- zoo(dat[, 4:6], order.by = dat$Date, frequency = 1)
  tsPQE
})
usethis::use_data(Canning, overwrite = TRUE)
