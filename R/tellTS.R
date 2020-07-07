# Call tell method for each column of ts.matrix, fun should return
#  long format data.frame: ts.id,variable,variables,...

tellTS <- function(x, ts.matrix, fun,
                   indices,
                   parallel = hydromad.getOption("parallel")[["tellTS"]],
                   ...) {
  UseMethod("tellTS")
}

tellTS.default <- function(x, ts.matrix, fun,
                           indices,
                           parallel = hydromad.getOption("parallel")[["tellTS"]],
                           ...) {
  force(fun)
  if (is.null(parallel)) parallel <- "none"
  if (missing(indices)) indices <- 1:ncol(ts.matrix)
  cat(sprintf(
    "Calculating sensitivity on %d data points %s parallelisation\n", length(indices),
    ifelse(parallel != "none", "*without*", "with")
  ))
  switch(parallel,
    "clusterApply" = {
      clusterEvalQ(cl, library(sensitivity))
      clusterEvalQ(cl, library(ff))
      clusterExport(cl, c("ts.matrix", "x", "fun"), envir = environment())
      results <- parLapply(cl, indices, function(i) {
        y <- ts.matrix[, i]
        tell(x, y, ...)
        return(fun(i, x))
      })
    },
    {
      results <- lapply(indices, function(i) {
        y <- ts.matrix[, i]
        tell(x, y, ...)
        return(fun(i, x))
      })
    }
  )
  return(do.call(rbind, results))
}

tellTS.sobol2002 <- function(x, ts.matrix, fun,
                             indices,
                             parallel = hydromad.getOption("parallel")[["tellTS"]],
                             ...) {
  if (missing(fun)) {
    fun <- function(i, x) {
      cbind(
        ts.id = i,
        variable = rownames(x$T),
        x$T[, c("original", "min. c.i.", "max. c.i.")]
      )
    }
  }

  tellTS.default(x, ts.matrix, fun,
    indices,
    parallel = hydromad.getOption("parallel")[["tellTS"]],
    ...
  )
}

tellTS.sobol2007 <- tellTS.sobol2002

tellTS.morris <- function(x, ts.matrix, fun,
                          indices,
                          parallel = hydromad.getOption("parallel")[["tellTS"]],
                          ...) {
  if (missing(fun)) {
    fun <- function(i, x) {
      res <- data.frame(
        ts.id = i,
        variable = colnames(x$ee),
        mu = apply(x$ee, 2, mean),
        mu.star = apply(x$ee, 2, function(x) mean(abs(x))),
        sigma = apply(x$ee, 2, sd)
      )
      rownames(res) <- NULL
      res
    }
  }

  tellTS.default(x, ts.matrix, fun,
    indices,
    parallel = hydromad.getOption("parallel")[["tellTS"]],
    ...
  )
}
