#' Rainfall and streamflow for Salmon Brook at Salmon Catchment.
#'
#' Daily rainfall and streamflow for Salmon Brook at Salmon Catchment (Western
#' Australia), from 1974-04-06 to 1999-03-16.  The catchment area is 0.83
#' square kilometers.
#'
#' Salmon Brook is located in the low-relief Darling Range of southwestern
#' Western Australia. Evapotranspiration consumes about 90 per cent of the
#' annual rainfall. The region is dominated by jarrah (Eucalyptus marginata)
#' forest. Surface soils are predominantly highly permeable sands and gravels
#' (Ye et al., 1997).
#'
#' \describe{ \item{Rainfall (P)}{ Daily rainfall (mm/day).  Rain gauge station
#' ID 509247 "SALMON BROOK @ SALMON CATCHMENT".  Latitude -33.4157; Longitude
#' 115.9846.  } \item{Streamflow (Q)}{ Daily mean streamflow (mm/day).  Stream
#' gauge ID 612011 "SALMON BROOK @ SALMON CATCHMENT".  Latitude -33.4176;
#' Longitude 115.9817.  } \item{Temperature (E)}{ Mean Maximum Temperature
#' Climate Data.  Product code: IDCJAC0002 reference: 00232595.
#'
#' Bureau of Meteorology station number: 9534 Station name: DONNYBROOK Latitude
#' -33.57; Longitude 115.82; Altitude 63m.  } }
#'
#' @name SalmonBrook
#' @docType data
#' @format A \code{\link{zoo}} object, of class \code{c("zooreg", "zoo")}.  It
#' is a regular time series indexed by days, in \code{Date} format.
#'
#' There are three columns, \code{P} (rainfall, mm/day) and \code{Q}
#' (streamflow, mm/day).  \code{E} (temperature, degrees C).
#' @seealso \code{\link{BinghamTrib}}, \code{\link{YeAl97}}
#' @references Ye, W., B. C. Bates, N. R. Viney, M. Sivapalan and A. J. Jakeman
#' (1997).  Performance of conceptual rainfall-runoff models in low-yielding
#' ephemeral catchments, \emph{Water Resources Research} 33, pp. 153--166.
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
#' data(SalmonBrook)
#' summary(SalmonBrook)
#' xyplot(SalmonBrook)
"SalmonBrook"
