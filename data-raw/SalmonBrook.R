SalmonBrook <- local({

  ## read files from doc directory
  pqdat <- read.table("SalmonBrook.dat",
    header = TRUE,
    col.names = c("Date", "P", "Q"), as.is = TRUE
  )
  tdat <- read.table("DonnybrookTemp.txt",
    header = TRUE,
    col.names = c("Year", "Month", "T")
  )
  ## convert dates
  pqdat$Date <- as.Date(pqdat$Date)
  tdat$Date <- with(tdat, as.Date(paste(Year, Month, 1, sep = "-")))

  ## convert from m^3/day to mm/day (i.e. divide by 1000 & catchment area km^2)
  pqdat$Q <- (pqdat$Q / 1000) / 0.83
  #    pqdat$Q <- convertFlow(pqdat$Q, from = "m^3", area.km2 = 0.83)

  ## zoo objects
  library(zoo)
  tsPQ <- zoo(pqdat[, 2:3], order.by = pqdat$Date, frequency = 1)
  tsT <- zoo(tdat$T, order.by = tdat$Date, frequency = 1)

  ## start temperature series from first month of PQ series
  tsT <- window(tsT, start = zoo:::as.Date.yearmon(as.yearmon(start(tsPQ))))
  ## make monthly temperature record a regular series (expand gaps)
  tsT <- merge(tsT, zoo(order.by = seq(start(tsT), end(tsT), by = "months")))
  ## fill monthly temperature gaps with seasonal averages
  months <- months(time(tsT))
  avgs <- tapply(coredata(tsT), months, mean, na.rm = TRUE)
  tsT[is.na(tsT)] <- avgs[match(months[is.na(tsT)], names(avgs))]

  tsPQE <- merge(tsPQ, E = tsT, all = TRUE)
  tsPQE$E <- na.locf(tsPQE$E, na.rm = FALSE)
  na.trim(tsPQE)
})
usethis::use_data(SalmonBrook, overwrite = TRUE)
