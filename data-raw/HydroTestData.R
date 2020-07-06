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
