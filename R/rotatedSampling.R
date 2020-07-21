## Sample from the rotated space defined by the eigenvectors of XtX


#' Sample within rotated feasible parameter space
#'
#' Using given feasible parameter set, rotate using eigenvectors and sample
#' within rotated hypercube
#'
#' Parameters may not be in feasible parameter set
#'
#' @param X matrix or data.frame of feasible parameters
#' @param samples Number of samples to take
#' @param expand Ratio with which to expand bounds in rotated parameter space,
#' as a buffer around the original hypercube
#' @param \dots Arguments to pass to \code{\link{parameterSets}}
#' @return data.frame of parameter values drawn from rotated hypercube
#' @author Joseph Guillaume
#' @seealso \code{\link{parameterSets}}
#' @references Karel, K. (1990). Membership-set estimation using random
#' scanning and principal component analysis. Mathematics and Computers in
#' Simulation 32(5-6): 535-543 DOI:
#' http://dx.doi.org/10.1016/0378-4754(90)90009-8.
#'
#' Also possibly related to Iman, R.L. and W.J. Conover (1982). A
#' distribution-free approach to inducing rank correlation among input
#' variables. Communications in Statistics - Simulation and Computation 11(3):
#' 311-334 DOI: http://dx.doi.org/10.1080/03610918208812265.
#' @examples
#'
#'
#' X <- matrix(c(1, 4, 2, 5, 1, 2, 2, 3), ncol = 2)
#' Y <- rotatedSampling(X, 1e3)
#' plot(X[, 1], X[, 2], col = "red", cex = 2)
#' points(Y[, 1], Y[, 2], pch = ".")
#' @export
rotatedSampling <- function(X, samples, expand = 0, ...) {
  if (is.null(colnames(X))) colnames(X) <- sprintf("X%d", 1:ncol(X))
  ## TODO: set fixed
  ## subtract mean for each parameter
  means <- colMeans(X)
  X2 <- t(t(X) - means)
  ## calculate dispersion matrix XtX square symmetrical
  disp <- t(X2) %*% X2
  ## eigenvalue decomposition
  e <- eigen(disp)
  V <- e$vectors
  ## rotate space using eigenvectors
  X3 <- X2 %*% V
  ## Calculate bounds of rotated space
  bounds <- apply(X3, 2, range)
  ## Optionally expand the bounds
  if (expand > 0) bounds <- apply(bounds, 2, function(x) c(x[1] * (1 - sign(x[1]) * expand), x[2] * (1 + sign(x[2]) * expand)))
  ## Sample from the bounds
  bounds <- lapply(apply(bounds, 2, list), unlist)
  names(bounds) <- colnames(X)
  Y <- parameterSets(bounds, samples, ...)
  ## Rotate back
  Ys <- t(t(as.matrix(Y) %*% t(V)) + means)
  Ys <- as.data.frame(Ys)
  names(Ys) <- names(Y)
  Ys
  ## TODO: adaptive, multiple rotations
}
