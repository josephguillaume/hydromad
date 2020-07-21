#' Rainfall and streamflow for Cotter River at Gingera.
#'
#' Daily rainfall and streamflow for Cotter River at Gingera (Australian
#' Capital Territory), from 1966-05-01 to 2003-06-07.  The catchment area is
#' 148 square kilometers.
#'
#' \describe{ \item{Rainfall (P)}{ Daily areal rainfall (mm/day).
#'
#' Derived from rain gauges operated by Bureau of Meteorology and EcoWise. An
#' area-weighted average was used, with weights determined from a long-term
#' spline-interpolated rainfall surface.  } \item{Streamflow (Q)}{ Daily mean
#' streamflow (mm/day).  Stream gauge ID 410730 "Cotter at Gingera".  Latitude
#' -33.592; Longitude 148.822.  } \item{Temperature (E)}{ Daily Maximum
#' Temperature.  Product code: IDCJAC0002 reference: 00232595.
#'
#' Bureau of Meteorology station number: 070014 Station name: CANBERRA AIRPORT
#' Latitude -35.305; Longitude 149.201; Altitude 578.4 m.  } }
#'
#' @name Cotter
#' @docType data
#' @format A \code{\link{zoo}} object, of class \code{c("zooreg", "zoo")}.  It
#' is a regular time series indexed by days, in \code{Date} format.
#'
#' There are three columns, \code{P} (areal rainfall, mm/day) and \code{Q}
#' (streamflow, mm/day).  \code{E} (temperature, degrees C).
#' @source EcoWise?
#'
#' Temperature: Copyright (c) Commonwealth of Australia.
#' @keywords datasets
#' @examples
#'
#' data(Cotter)
#' summary(Cotter)
#' xyplot(Cotter)
"Cotter"
