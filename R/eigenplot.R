eigen.plot.single <- function(e, max.value = NA) {
  stopifnot(!is.null(e$values))
  stopifnot(!is.null(e$vectors))
  evs <- sqrt(abs(e$values))
  evecs <- e$vectors
  a <- evs[1]
  b <- evs[2]
  x0 <- 0
  y0 <- 0
  alpha <- atan(evecs[, 1][2] / evecs[, 1][1])
  theta <- seq(0, 2 * pi, length = (1000))
  x <- x0 + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
  y <- y0 + a * cos(theta) * sin(alpha) + b * sin(theta) * cos(alpha)
  ## plot(x, y, type = "l", asp = 1, ann = FALSE, axes=FALSE, xlim = c(-0.14, 0.14), ylim = c(-0.14, 0.14))
  if (!is.na(max.value)) {
    plot(x, y, type = "l", asp = 1, ann = FALSE, axes = FALSE, xlim = c(-max.value, max.value), ylim = c(-max.value, max.value))
  } else {
    plot(x, y, type = "l", asp = 1, ann = FALSE, axes = FALSE)
  }
  arrows(x0 = 0, y0 = 0, x1 = a * evecs[, 1][1], y1 = a * evecs[, 1][2], length = 0)
  arrows(x0 = 0, y0 = 0, x1 = b * evecs[, 2][1], y1 = b * evecs[, 2][2], length = 0)
}

eigen.plot <- function(obj, fixed.axis = TRUE) {
  stopifnot(inherits(obj, "rsm"))
  ## Calculated pair-wise eigen values first to allow fixed axis
  Bm <- list() # B matrix of coefficients
  inv.Bm <- list()
  eg <- list() # eigen values and vectors
  npar <- nrow(obj$B)
  for (i in 1:(npar - 1)) {
    for (j in 2:npar) {
      if (i >= j) next
      Bm[[paste(i, j)]] <- obj$B[c(i, j), c(i, j)]
      inv.Bm[[paste(i, j)]] <- -(solve(Bm[[paste(i, j)]]))
      eg[[paste(i, j)]] <- eigen(inv.Bm[[paste(i, j)]])
    }
  }

  fixed.axis <- ifelse(fixed.axis,
    sqrt(max(sapply(eg, function(x) max(abs(x$values))))) * 1.1,
    NA
  )
  par(mfcol = c(npar - 1, npar - 1))
  par(oma = c(0, 4, 4, 0)) ## for text for axis
  par(bty = "n") ## remove border
  par(mar = c(0, 0, 0, 0)) ## remove useless margin
  for (i in 1:(npar - 1)) {
    for (j in 2:npar) {
      if (i >= j) {
        frame()
        if (j == 2) mtext(colnames(obj$B)[i])
        next
      }
      ## eigen.plot(eigen(-solve(obj$B[c(i,j),c(i,j)]))) ##if not precomputed
      eigen.plot.single(eg[[paste(i, j)]], fixed.axis)
      if (j == 2) mtext(colnames(obj$B)[i]) ## top label
      if (i == 1) mtext(colnames(obj$B)[j], side = 2) ## side label
    }
  }
  invisible(eg)
}
