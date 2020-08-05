## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Joseph Guillaume <josephguillaume@gmail.com>
##



#' Evaluate a model for a matrix of parameters
#'
#' \code{evalPars} evaluates a model for a named matrix of parameters.
#'
#' \code{evalPars} is used in conjunction with \code{getFreeParsRanges}
#' and the \code{sensitivity} package to perform sensitivity analysis. See
#' \code{demo(sensitivity)}. Note that the objective function may more
#' generally return any scalar result, e.g. a scalar prediction calculated only
#' using X.
#'
#' Individual model evaluations are generally very fast, so parallelisation is
#' only really worthwhile when large numbers of evaluations are needed.
#' \code{"clusterApply"} has a slightly lower overhead. \code{"foreach"} allows
#' a broader range of options.
#'
#' @name evalPars
#' @aliases getFreeParsRanges
#' @param par.matrix Named matrix or data.frame of parameter values, with each
#' row corresponding to a model realisation to evaluate
#' @param object an object of class \code{hydromad}.
#' @param objective the objective function or expression, which can refer to Q
#' and X.  See \code{\link{objFunVal.hydromad}}
#' @param parallel Which parallelisation method to use. Options are
#' \code{"clusterApply"} and \code{"foreach"}. For any other value, parameters
#' will be evaluated sequentially. Default is \code{"none"}. For more
#' information see \code{\link{hydromad_parallelisation}}.
#' @return evalPars returns a vector of objective function values,
#' corresponding to the evaluation of each row of \code{par.matrix}.
#' @author Joseph Guillaume
#' @seealso \code{\link{objFunVal.hydromad}}, \code{\link{parameterSets}},
#' \code{getFreeParsRanges}
#' @keywords models
#' @examples
#'
#' data(Cotter)
#' obs <- Cotter[1:1000]
#' modx <- hydromad(obs,
#'   sma = "cmd", routing = "expuh",
#'   tau_s = c(2, 100), v_s = c(0, 1)
#' )
#'
#' ## Sample 10 random parameter sets for parameters with defined ranges
#' pars <- parameterSets(getFreeParsRanges(modx), 10, method = "random")
#'
#' ## Return the default objective function value for each model realisation
#' evalPars(pars, object = modx)
#'
#' ## Calculate the 20%ile flow for each model realisation
#' evalPars(pars, object = modx, objective = ~ quantile(X, 0.2))
#'
#' ## Alternatively, sample 10 random parameter sets from all parameters
#' ##  This allows specifying discrete values of parameters
#' pars <- parameterSets(coef(modx, warn = FALSE), 10, method = "random")
#' @export
evalPars <- function(par.matrix, object, objective = hydromad.getOption("objective"),
                     parallel = hydromad.getOption("parallel")[["evalPars"]]) {
  stopifnot(inherits(object, "hydromad"))

  # Sets default settings for parallelisation if missing
  parallel <- hydromad.parallel(parallel)

  switch(parallel$method,
    "foreach" = {
      opts <- hydromad.options()
      export <- parallel$export
      objs <- foreach::foreach(
        p = iterators::iter(par.matrix, by = "row"),
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


#' @export
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

#' @import utils
utils::globalVariables(c("cl", "p"))
