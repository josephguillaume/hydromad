#' Rainfall and streamflow for Murrindindi River at Murrindindi above Colwells.
#'
#' Daily rainfall and streamflow for Murrindindi River at Murrindindi above
#' Colwells (Victoria, Australia), from 1975-06-08 to 1998-06-30.  The
#' catchment area is 105.9 square kilometers.
#'
#' \describe{ \item{Rainfall (P)}{ Daily areal rainfall (mm/day).
#'
#' Unknown...
#'
#' Derived from rain gauges operated by Bureau of Meteorology.  An
#' area-weighted average was used, with weights determined from a long-term
#' spline-interpolated rainfall surface.  } \item{Streamflow (Q)}{ Daily mean
#' streamflow (mm/day).  Stream gauge ID 405205A "Murrindindi @ Murrindidi
#' above Colwells".  Latitude -37.412; Longitude 145.558.  } \item{Temperature
#' (E)}{ Unknown.  } }
#'
#' @name Murrindindi
#' @docType data
#' @format A \code{\link{zoo}} object, of class \code{c("zooreg", "zoo")}.  It
#' is a regular time series indexed by days, in \code{Date} format.
#'
#' There are three columns, \code{P} (areal rainfall, mm/day) and \code{Q}
#' (streamflow, mm/day).  \code{E} (temperature, degrees C).
#' @source Department of Sustainability and Environment (VIC).
#'
#' Temperature: Copyright (c) Commonwealth of Australia.
#' @keywords datasets
#' @examples
#'
#' data(Murrindindi)
#' summary(Murrindindi)
#' xyplot(Murrindindi)
"Murrindindi"
