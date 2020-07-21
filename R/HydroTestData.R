#' A simple simulated dataset for use in testing hydrological models.
#'
#' A simple simulated dataset for use in testing hydrological models.
#'
#' \describe{ \item{Rainfall (P)}{ a regular series of impulses (every 20 time
#' steps). Each of these pulses have a value of 5, except one which has a value
#' of 20.  } \item{Temperature (E)}{ a sine wave ranging from 0 to 30.  }
#' \item{Streamflow (Q)}{ proportional to the square of rainfall and inversely
#' to temperature, then filtered with a second-order autoregressive
#' \code{\link{filter}}.  } }
#'
#' @format A \code{\link{zooreg}} object with 730 (365 * 2) time steps.
#'
#' There are three columns, \code{P} (simulated areal rainfall, mm/day),
#' \code{E} (simulated temperature, degrees Celcius) and \code{Q} (simulated
#' streamflow, mm/day).
#' @keywords datasets
#' @examples
#'
#' data(HydroTestData)
#' summary(HydroTestData)
#' xyplot(HydroTestData)
"HydroTestData"
