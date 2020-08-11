#' @md
#' @name Corin
#' @title Daily dataset for the Corin catchment.
#' @description A daily dataset containing precipitation, streamflow, potential
#' evapotranspiration and average temperature for the Corin catchment (147 km²),
#' ACT, Australia. The Corin catchment is the headwaters of the greater Cotter
#' River catchment, which is a tributary to the Murrumbidgee River.
#' @usage data("Corin")
#' @format
#' A `zoo` object of class `"zooreg", "zoo"`. It is a regular time series
#' indexed by days, in `Date` format.
#' There are four columns:
#' * `P` (precipitation, mm/day)
#' * `Q` (streamflow, mm/day)
#' * `E` (potential evapotranspiration, mm/day)
#' * `T` (daily average temperature, ºC)
#' @details
#' * Daily rainfall at Cotter Hut (ID: 570946).
#' * Daily mean streamflow at Cotter River at Gingera (ID: 410730).
#' * Daily Morton's potential ET calculated using Corin Dam (ID: 570947) weather
#' station data.
#' * Daily average temperature as measured at Corin Dam weather station
#' (ID: 570947).
#' @source
#' Data obtained through the Australian Bureau of Meteorology Water Data web services API.
#' Data owner: ACT - Icon Water Limited.
#' @examples
#' data(Corin)
#' summary(Corin)
#' xyplot(Corin)
#' @keywords datasets
"Corin"
