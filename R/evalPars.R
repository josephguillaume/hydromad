## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Joseph Guillaume <josephguillaume@gmail.com>
##

evalPars <- function(par.matrix, object, objective = hydromad.getOption("objective"),
                     parallel = hydromad.getOption("parallel")[["evalPars"]]) {
  stopifnot(inherits(object, "hydromad"))

  # Sets default settings for parallelisation if missing
  parallel <- hydromad.parallel(parallel)

  switch(parallel$method,
    "foreach" = {
      opts <- hydromad.options()
      export <- parallel$export
      objs <- foreach(
        p = iter(par.matrix, by = "row"),
        .packages = parallel$packages,
        .inorder = TRUE,
        .export = parallel$export,
        .final = function(x) as.numeric(x),
        .options.redis = list(async = parallel$async)
      ) %dopar% {
        # Work-around for hydromad functions to have access to .export
        for (e in export) assign(e, get(e), envir = .GlobalEnv)
        # Work-around to use same opts as in user's environment
        hydromad.options(opts)
        thisMod <- update(object, newpars = p)
        if (!isValidModel(thisMod)) {
          return(NA)
        }
        objFunVal(thisMod, objective = objective)
      }
    },
    "clusterApply" = {
      if (length(parallel$packages) > 0) lapply(parallel$packages, function(pkg) clusterCall(cl, library, pkg, character.only = TRUE))
      if (length(parallel$export) > 0) clusterExport(cl, parallel$export)
      objs <- parApply(
        cl = cl, par.matrix, 1,
        function(p, object, objective) {
          thisMod <- update(object, newpars = p)
          if (!isValidModel(thisMod)) {
            return(NA)
          }
          objFunVal(thisMod, objective)
        }, object = object, objective = objective
      )
    },
    objs <- apply(par.matrix, 1, function(p) {
      thisMod <- update(object, newpars = p)
      if (!isValidModel(thisMod)) {
        return(NA)
      }
      objFunVal(thisMod, objective)
    })
  ) ## switch parallel
  return(objs)
}


getFreeParsRanges <- function(object) {
  stopifnot(inherits(object, "hydromad"))
  ## identify varying parameters
  par.ranges <- suppressWarnings(coef(object))
  free <- sapply(par.ranges, function(x) {
    !inherits(x, "AsIs") && length(x) == 2 && (diff(range(x)) >
      0)
  })
  par.ranges[free]
}
