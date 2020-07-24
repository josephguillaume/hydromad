## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##




#' Estimate the dead time between input and output
#'
#' Use cross-correlation to estimate the delay between an input time series and
#' (rises in) the corresponding output time series.
#'
#' The estimated delay is the one maximising the cross-correlation function.
#'
#' @importFrom stats na.exclude ccf frequency
#'
#' @param DATA a \code{\link{ts}}-like object with named components: \describe{
#' \item{list("U")}{ input (forcing) time series. } \item{list("Q")}{ output
#' (response) time series. } }
#' @param rises use only rises in the output to estimate delay.
#' @param lag.max largest delay (in time steps) to consider.
#' @param n.estimates number of estimates of delay to produce.
#' @param negative.ok to allow negative delay times to be considered.
#' @param na.action handler for missing values.  The default removes leading
#' and trailing NAs only.  Use \code{na.exclude} to remove all NAs, but the
#' result will not be a valid autocorrelation sequence.
#' @param plot plot the cross-correlation function.
#' @param main title for plot.
#' @param \dots further arguments passed to \code{\link{ccf}} or on to
#' \code{\link{plot.acf}}.
#' @return The estimated delay as an integer number of time steps.  If
#' \code{n.estimates > 1}, that number of integer delays, ordered by the CCF.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{ccf}},\code{\link{estimateDelayFrac}}
#' @keywords ts
#' @examples
#'
#' set.seed(1)
#' x <- ts(pmax(0, rgamma(200, shape = 0.1, scale = 20) - 5))
#' ## simulate error as multiplicative uniform random
#' y <- x * runif(200, min = 0.5, max = 1.5)
#' ## and resample 10 percent of time steps
#' ii <- sample(seq_along(y), 20)
#' y[ii] <- rev(y[ii])
#' ## apply recursive filter and lag
#' y <- filter(y, 0.8, method = "r")
#' y <- lag(y, -2) # true delay is 2
#' plot(ts.union(y, x))
#' ## based on cross correlation function:
#' estimateDelay(ts.union(y, x), rises = FALSE, plot = TRUE)
#' ## based on ccf with flow rises only:
#' estimateDelay(ts.union(y, x), plot = TRUE)
#' @export
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
