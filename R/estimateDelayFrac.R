#' Estimate the dead time between input and output, with a fractional component
#' (redistribution of the input)
#'
#' Optimises the delay \code{TDopt} using \code{\link{lagFrac}} to maximise the
#' correlation between the delayed input time series and (rises in) the
#' corresponding time series
#'
#' \code{\link{estimateDelay}} may be preferred if there's no good reason to
#' have a fractional lag/redistribution of the input across days.
#'
#' @importFrom stats optimise
#'
#' @param DATA a \code{\link{ts}}-like object with named components: \describe{
#' \item{list("U")}{ input (forcing) time series. } \item{list("Q")}{ output
#' (response) time series. } }
#' @param rises use only rises in the output to estimate delay.
#' @param lag.max largest delay (in time steps) to consider.
#' @return The estimated delay as an integer number of time steps.
#' @author Joseph Guillaume
#' @seealso \code{\link{estimateDelay}}, \code{\link{lagFrac}}
#' @keywords ts
#' @examples
#'
#' L <- 0.6 ## Lag of 0.6
#' P <- c(2, 0, 5, 1, 6, 10, 0, 0, 0)
#' V1 <- lagFrac(P, L)
#'
#' estimateDelay(cbind(P, V1), rises = FALSE)
#' estimateDelayFrac(cbind(U = P, Q = V1), lag.max = 5, rises = FALSE)
#' @export
estimateDelayFrac <- function(DATA, rises = TRUE, lag.max = hydromad.getOption("max.delay")) {
  DATA <- as.ts(DATA)
  if (NROW(DATA) <= 1) {
    return(NA_integer_)
  }
  iQ <- 1
  if ("Q" %in% colnames(DATA)) {
    iQ <- which("Q" == colnames(DATA))[1]
  }
  Q <- DATA[, iQ]
  U <- DATA[, if (iQ == 1) {
    2
  } else {
    1
  }]
  if (all(Q == 0, na.rm = TRUE)) {
    return(NA_integer_)
  }
  do.rises <- rises
  if (do.rises) {
    rises <- function(x) {
      x[] <- c(0, pmax(diff(x), 0))
      x
    }
    Q <- rises(Q)
  }
  return(optimise(function(L) cor(lagFrac(U, L), Q, use = "complete.obs"), interval = c(0, lag.max), maximum = T)$maximum)
}

#################################################################################
# x<-mod$data
# event <- eventseq(x$P, thresh = 10, inthresh = 1, indur = 7, continue = TRUE)
# event.1lag <- lag(event, 1) # to consider the increment of the first obsQ
# whole.lag<-estimateDelayFrac(x)
#
# PQEM=mod$data
# lag.rain<-eventapply(PQEM,event.1lag,by.column=FALSE,function(x){
# ans.lag<-estimateDelayFrac(x, lag.max = 4)
# if (is.na(ans.lag)) {ans.lag<-whole.lag} # if obsQ=0 then delay=NA, put whole delay value
# x2<-merge(x, P1 = lagFrac(x[,1], ans.lag), ans.lag, all = c(TRUE, FALSE))
# return(x2)
# })

# require(reshape)
# melt.lag<-melt(lag.rain)
# write.csv(melt.lag, file = "lagrain.csv")
