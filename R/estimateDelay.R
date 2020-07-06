## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


estimateDelay <-
  function(DATA = data.frame(U = , Q = ),
           rises = TRUE,
           lag.max = hydromad.getOption("max.delay"),
           n.estimates = 1, negative.ok = FALSE,
           na.action = na.exclude, plot = FALSE, main = NULL,
           ...) {
    ## get data into the right form
    DATA <- as.ts(DATA)
    DATA <- na.action(DATA)
    if (NROW(DATA) <= 1) {
      return(NA_integer_)
    }
    ## which column represents streamflow?
    iQ <- 1
    if ("Q" %in% colnames(DATA)) {
      iQ <- which("Q" == colnames(DATA))[1]
    }
    Q <- DATA[, iQ]
    U <- DATA[, if (iQ == 1) 2 else 1]

    if (all(Q == 0, na.rm = TRUE)) {
      return(NA_integer_)
    }

    ylab <- "CCF"
    do.rises <- rises
    if (do.rises) {
      ## backwards difference, i.e. rises are placed at their end time
      rises <- function(x) {
        x[] <- c(0, pmax(diff(x), 0))
        x
      }
      Q <- rises(Q)
    }
    if (is.null(main)) {
      main <- "Cross-correlation"
      if (do.rises) main <- paste(main, "with rises only")
    }
    if (sum(complete.cases(Q, U)) == 0) {
      return(NA_integer_)
    }
    ans <- ccf(Q, U,
      lag.max = lag.max, na.action = na.action,
      plot = plot, main = main, ...
    )
    ## ans$lag[which.max(ans$acf)]
    est <- ans$lag[order(ans$acf, decreasing = TRUE)]
    ## omit any negative delays or NAs
    invalid <- is.na(est)
    if (negative.ok == FALSE) {
      invalid <- invalid | (est < 0)
    }
    est <- est[!invalid]
    if (length(est) == 0) {
      return(NA_integer_)
    }
    est <- est[seq(n.estimates)]
    ## convert from "basic units" to time steps
    est <- est * frequency(Q)
    est
  }
