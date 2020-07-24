

## Time series ('ts' or 'zoo') comparison on corresponding times
ts_equals <- function(expected, actual, start = NULL, end = NULL, trim = FALSE) {
  function(actual) {
    windowts <-
      window(ts.intersect(as.ts(expected), as.ts(actual)),
        start = start, end = end
      )
    if (trim) {
      windowts <- na.trim(windowts)
    }
    identical(windowts[, 1])(windowts[, 2])
  }
}


## Time series ('ts' or 'zoo') comparison on corresponding times
ts_equals2 <- function(expected, actual, start = NULL, end = NULL, trim = FALSE) {
  windowts <-
    window(ts.intersect(as.ts(expected), as.ts(actual)),
      start = start, end = end
    )
  if (trim) {
    windowts <- na.trim(windowts)
  }
  windowts[, 1]
}
