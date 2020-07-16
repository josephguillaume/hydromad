#' Rainfall and streamflow for Bingham River Trib at Ernies Catchment.
#' 
#' Daily rainfall and streamflow for Bingham River Trib at Ernies Catchment
#' (Western Australia), from 1974-05-18 to 2008-11-02.  The catchment area is
#' 2.68 square kilometers.
#' 
#' Salmon Brook is located in the low-relief Darling Range of southwestern
#' Western Australia.
#' 
#' \describe{ \item{Rainfall (P)}{ Daily rainfall (mm/day).  Rain gauge station
#' ID 509249 "BINGHAM RIVER TRIB @ ERNIES CATCHMENT".  Latitude -33.2921;
#' Longitude 116.4451.  } \item{Streamflow (Q)}{ Daily mean streamflow
#' (mm/day).  Stream gauge ID 612008 "BINGHAM RIVER TRIB @ ERNIES CATCHMENT".
#' Latitude -33.2939; Longitude 116.4449.  } \item{Temperature (E)}{ Mean
#' Maximum Temperature Climate Data.  Product code: IDCJAC0002 reference:
#' 00232595.
#' 
#' Bureau of Meteorology station number: 9534 Station name: DONNYBROOK Latitude
#' -33.57; Longitude 115.82; Altitude 63m.  } }
#' 
#' @name BinghamTrib
#' @docType data
#' @format A \code{\link{zoo}} object, of class \code{c("zooreg", "zoo")}.  It
#' is a regular time series indexed by days, in \code{Date} format.
#' 
#' There are three columns, \code{P} (rainfall, mm/day) and \code{Q}
#' (streamflow, mm/day).  \code{E} (temperature, degrees C).
#' @seealso \code{\link{SalmonBrook}}
#' @source Water INformation (WIN) database - discrete sample data.
#' [2009-07-09]. Department of Water, Water Information Provision section,
#' Perth Western Australia.
#' 
#' Hydstra database - time-series data.  [2009-07-09].  Department of Water,
#' Water Information Provision section, Perth Western Australia.
#' 
#' Copyright (c) State of Western Australia (Department of Water).
#' 
#' Temperature: Copyright (c) Commonwealth of Australia. Created on Tue 07 Jul
#' 2009 05:48:38 AM GMT from Climate Data Online, Bureau of Meteorology.
#' http://www.bom.gov.au/climate/averages
#' @keywords datasets
#' @examples
#' 
#' data(BinghamTrib)
#' summary(BinghamTrib)
#' xyplot(BinghamTrib)
#' 
BinghamTrib <- local({
  
  ## read files from doc directory
  pqdat <- read.table("data-raw/BinghamTrib.dat",
                      header = TRUE,
                      col.names = c("Date", "P", "Q"), as.is = TRUE
  )
  tdat <- read.table("data-raw/DonnybrookTemp.txt",
                     header = TRUE,
                     col.names = c("Year", "Month", "T")
  )
  ## convert dates
  pqdat$Date <- as.Date(pqdat$Date)
  tdat$Date <- with(tdat, as.Date(paste(Year, Month, 1, sep = "-")))
  
  ## convert from m^3/day to mm/day (i.e. divide by 1000 & catchment area km^2)
  pqdat$Q <- (pqdat$Q / 1000) / 2.68
  #    pqdat$Q <- convertFlow(pqdat$Q, from = "m^3", area.km2 = 2.68)
  
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
usethis::use_data(BinghamTrib, overwrite = TRUE)
