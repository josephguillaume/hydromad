#' Non-integer time delay
#'
#' Delays P by a non-integer number of timesteps using a filter
#'
#'
#' @param P Input to delay (typically rainfall)
#' @param TDopt Positive number with which to delay \code{P}
#' @return vector of same length as P
#' @author Thanks to Charles Perrin for original code.  Implemented in R by
#' Joseph Guillaume
#' @seealso \code{\link{estimateDelayFrac}}
#' @keywords ts
#' @examples
#'
#' L <- 0.6
#' P <- c(2, 0, 5, 1, 6, 10, 0, 0, 0)
#' V1 <- lagFrac(P, L)
#' @export
lagFrac <- function(P, TDopt) {
  TDopt <- TDopt + 1
  I <- floor(TDopt)
  F <- TDopt - I
  X <- rep(0, ceiling(TDopt))
  X[I] <- 1 - F
  X[I + 1] <- F
  filter.pad0 <- function(x, f) {
    y <- x
    y[] <- filter(c(rep(0, length(f)), x),
      filter = f, sides = 1
    )[-(1:length(f))]
    y
  }
  V <- filter.pad0(P, X)
  if (I > 1) V[1:(I - 1)] <- ifelse(V[1:(I - 1)] == 0, NA, V[1:(I - 1)])
  return(V)
}
