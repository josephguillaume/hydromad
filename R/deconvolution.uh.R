#' Unit Hydrograph using deconvolution
#'
#' Estimates the unit hydrograph using deconvolution.
#'
#' Estimates the unit hydrograph by fourier transform deconvolution of the
#' ratio of the cross correlation and auto-correlation functions.
#'
#' @importFrom stats ccf na.pass acf dnorm fft
#' @importFrom graphics plot
#'
#' @param P Rainfall or effective rainfall time series
#' @param Q Flow time series
#' @param FWHM Full Width Half Maximum to use for Gaussian apodisation function
#' @param do.plot Whether to plot the unit hydrograph immediately
#' @return Vector of length \code{(length(P)-1)*2} representing deconvolved
#' unit hydrograph. Note that entries from \code{length(P)} represent are
#' obtained by negative lags.  Optionally plot the unit hydrograph.
#' @author Joseph Guillaume
#' @seealso \code{\link{expuh},\link{lambda},\link{armax},\link{powuh}}
#' @references Barry
#' @keywords ts
#' @examples
#'
#' data(Murrindindi)
#' h <- deconvolution.uh(Murrindindi$P, Murrindindi$Q, do.plot = TRUE)
#' head(h)
#' @export
deconvolution.uh <- function(P, Q, FWHM = length(P), do.plot = FALSE) {
  ## Calculate PQ and PP cross correlations
  c <- ccf(Q, P, lag.max = length(P), plot = FALSE, na.action = na.pass)
  c <- c$acf[, 1, 1]
  a <- acf(P, lag.max = length(P), plot = FALSE, na.action = na.pass)
  a <- c(rev(a$acf[, 1, 1]), a$acf[-1, 1, 1])

  ## Apply apodisation function
  apod <- dnorm(((-length(P) + 1):(length(P) - 1)), sd = FWHM / 2 / sqrt(2 * log(2)))
  c <- c * apod
  a <- a * apod

  ## Deconvolve Q=P*H
  H <- fft(c) / fft(a)
  h <- fft(H, inverse = T) / length(c)
  stopifnot(Im(h) < 1e-10)
  h <- Re(h)
  if (do.plot) {
    n0 <- 10
    n1 <- 100
    plot(c(-n0:-1, 0:n1), c(h[(length(h) - n0 + 1):length(h)], h[1:(n1 + 1)]),
      col = "green", type = "l",
      xlab = "day", ylab = "Q", ylim = c(-0.1, 1)
    )
  }
  invisible(h)
}
