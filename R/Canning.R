#' Rainfall, streamflow and potential evaporation data for Canning River at
#' Scenic Drive.
#' 
#' Daily rainfall, streamflow and potential evaporation for Canning River at
#' Scenic Drive upstream of Canning Dam (Western Australia), from 1977-01-01 to
#' 1987-12-31.  The catchment area is 517 square kilometers.
#' 
#' The Canning River is a tributary of the Swan River, and is located South
#' East of Perth, Western Australia.
#' 
#' \describe{ \item{Rainfall (P)}{ Daily areal rainfall (mm/day).  Origin
#' unknown, probably spatial interpolation by Barry Croke.  } \item{Streamflow
#' (Q)}{ Daily mean streamflow (mm/day).  Stream gauge ID 616024 "Canning River
#' at Scenic Drive".  Latitude -33.4176; Longitude 115.9817.  } \item{Potential
#' Evapotranspiration (E)}{ Origin Unknown.} }
#' 
#' @name Canning
#' @docType data
#' @format A \code{zoo} object, of class \code{"zooreg", "zoo"}.  It
#' is a regular time series indexed by days, in \code{Date} format.
#' There are three columns, \code{P} (areal rainfall, mm/day) and \code{Q}
#' (streamflow, mm/day). \code{E} (potential evapotranspiration, mm/day).
#' 
#' @source Talk to Barry Croke...
#' @keywords datasets
#' @examples
#' 
#' data(Canning)
#' summary(Canning)
#' xyplot(Canning)
#' 
Canning <- local({
  dat <- read.table("data-raw/Canning.csv",
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
