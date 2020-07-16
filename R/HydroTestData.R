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
#' 
HydroTestData <- local({
  timeseq <- seq(as.POSIXct("2000-01-01", tz = "GMT"),
                 as.POSIXct("2000-03-31", tz = "GMT"),
                 by = "3 hours"
  )
  len <- length(timeseq)
  ## regular rainfall impulse
  P <- ts(rep(0, length = len))
  P[seq(10, len - 1, by = 24)] <- 6
  ## make one rain event much larger
  P[quantile(which(P > 0), 0.67, type = 1)] <- 24
  ## sine wave for temperature
  E <- ts(15 + 15 * sin(seq(0, 4 * pi, length = len)))
  ## flow based on square of rainfall and inverse to temperature
  Q <- stats::filter(0.01 * P^2 * (1 - E / max(E)),
                     filter = c(1.4, -0.45), method = "recursive"
  )
  as.zooreg(zoo(cbind(P = P, E = E, Q = Q), order.by = timeseq))
})
usethis::use_data(HydroTestData, overwrite = TRUE)
