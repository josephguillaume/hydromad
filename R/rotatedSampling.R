## Sample from the rotated space defined by the eigenvectors of XtX
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
