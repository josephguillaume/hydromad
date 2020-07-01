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
