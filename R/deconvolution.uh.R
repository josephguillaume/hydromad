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
