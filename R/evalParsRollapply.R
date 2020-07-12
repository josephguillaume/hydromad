## Run model and return a vector (usually a timeseries)
## Scalable in length of ts and number of parameters
## Parallelised with results stored on disk
## tempfile is automatically deleted when exiting from R
evalParsTS <- function(par.matrix, object,
                       fun = function(thisMod) fitted(thisMod),
                       length.out = NULL,
                       ...,
                       parallel = hydromad.getOption("parallel")[["evalParsTS"]],
                       filehash.name = tempfile()) {
  .args <- list(...)

  if (is.null(length.out)) {
    warning("length.out missing, running model to evaluate it")
    res <- fun(update(object, newpars = par.matrix[1, ]), ...)
    length.out <- length(res)
    cat("fun returns vector with length.out=", length.out, "\n")
  }

  # Sets default settings for parallelisation if missing
  parallel <- hydromad.parallel(parallel)

  if (parallel$method == "foreach" && !is.null(filehash.name)) {
    filehash.name <- NULL
    warning("ignoring filehash.name, 'foreach' parallelisation does not support writing directly to disk")
  }

  if (is.null(filehash.name)) {
    ## Don't use disk (ff)
    results <- matrix(NA, nrow = nrow(par.matrix), ncol = length.out)
    if (parallel$method == "clusterApply") {
      warning("setting parallel$method='none', 'clusterApply' parallelisation requires filehash.name to be non-null")
      parallel$method <- "none"
    }
  } else {
    ## Use disk
    if (!require("parallel")) stop("package parallel is required for evalParsTS if filehash.name is not NULL and parallel$method is not 'foreach'")
    if (!require("ff")) stop("package ff is required for evalParsTS if filehash.name is not NULL and parallel$method is not 'foreach'")
    results <- ff::ff(vmode = "double", dim = c(nrow(par.matrix), length.out), filename = filehash.name)
  }

  ## Do runs, storing all ts
  cat(sprintf("Running %d model evaluations with parallelisation='%s'\n", nrow(par.matrix), parallel$method))
  switch(parallel$method,
    "foreach" = {
      opts <- hydromad.options()
      export <- parallel$export
      results <- foreach::foreach(
        p = iter(par.matrix, by = "row"),
        .packages = parallel$packages,
        .inorder = TRUE,
        .export = export,
        .final = function(x) do.call(rbind, lapply(x, coredata)),
        .options.redis = list(async = parallel$async)
      ) %dopar% {
        # Work-around for hydromad functions to have access to .export
        for (e in export) assign(e, get(e), envir = .GlobalEnv)
        # Work-around to use same opts as in user's environment
        hydromad.options(opts)
        thisMod <- update(object, newpars = p)
        if (!isValidModel(thisMod)) {
          return(NULL)
        }
        do.call(fun, modifyList(.args, list(thisMod = thisMod)))
      }
    },
    "clusterApply" = {
      lapply(c("ff", parallel$packages), function(pkg) parallel::clusterCall(cl, library, pkg, character.only = TRUE))
      if (length(parallel$export) > 0) parallel::clusterExport(cl, parallel$export)
      parallel::clusterExport(cl, c("object"), envir = environment())
      parallel::parLapply(
        cl = cl, as.list(1:nrow(par.matrix)),
        function(ip) {
          thisMod <- update(object, newpars = par.matrix[ip, ])
          if (!isValidModel(thisMod)) {
            return(NULL)
          }
          results[ip, ] <- do.call(fun, modifyList(.args, list(thisMod = thisMod)))
          ip
        }
      )
    },
    {
      lapply(
        as.list(1:nrow(par.matrix)),
        function(ip) {
          thisMod <- update(object, newpars = par.matrix[ip, ])
          if (!isValidModel(thisMod)) {
            return(NULL)
          }
          results[ip, ] <<- do.call(fun, modifyList(.args, list(thisMod = thisMod)))
          ip
        }
      )
    }
  ) ## switch parallel
  ## ff matrix of ts
  results
}

## Calculate objective function on a rolling window
evalParsRollapply <- function(par.matrix, object,
                              width = 30,
                              objective = hydromad.getOption("objective"),
                              parallel = hydromad.getOption("parallel")[["evalParsTS"]],
                              filehash.name = tempfile()) {
  fun <- function(thisMod, width, objective) {
    rollapply(cbind(Q = observed(thisMod), X = fitted(thisMod)),
      width = width, by.column = FALSE,
      FUN = objFunVal, objective = objective
    )
  }

  evalParsTS(par.matrix, object, fun,
    length.out = length(observed(object)) - width + 1,
    parallel = parallel,
    filehash.name = filehash.name,
    width = width, objective = objective
  )
}
