## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##




#' Calculate objective function value for a fitted model.
#'
#' Calculate objective function value for a fitted model.
#'
#'
#' The objective function is given as a function with arguments \code{Q},
#' \code{X} and \code{\dots{}}, and optionally other arguments.  \code{Q} and
#' \code{X} represent observed and modelled flow, respectively.  It should
#' return a single numeric value.
#'
#' For more advanced use it may also refer to arguments \code{U} (modelled
#' effective rainfall), or \code{P} (observed rainfall), and more generally it
#' may refer to \code{model}, and so may extract other items of data,
#' parameters, etc.
#'
#' The default (unless changed in \code{hydromad.options("objective")}) is a
#' weighted sum of the R Squared (coefficient of determination) of square-root
#' transformed data, and the relative bias.
#'
#' See \code{\link{hydromad.stats}} for examples of how to specify other fit
#' statistics.
#'
#' @importFrom parallel parSapply
#'
#' @aliases objFunVal objFunVal.hydromad objFunVal.default objFunVal.runlist
#' @param x object from which to calculate stats. For the \code{hydromad}
#' method, this should be a fitted \code{hydromad} model object, i.e. it must
#' have specific parameter values, not ranges. For the default method, this
#' should be a matrix-like object with named columns \code{Q}, \code{X} and
#' optionally \code{U} and \code{P}.
#' @param objective the objective function, or a list of objective functions.
#' See Details.
#'
#' For a \code{runlist}, either: 1) a function taking a runlist as first
#' argument or 2) a list with two elements, the first as for \code{hydromad}
#' and \code{default}, the second a function that aggregates the objectives to
#' a single value.
#' @param ...  ignored for \code{hydromad} and \code{default}, passed to
#' \code{objFunVal.hydromad} for \code{runlist}
#' @param all passed to \code{fitted} and \code{observed}.
#' @param nan.ok by default, an error is thrown if the result is \code{NaN}.
#' Set this argument to \code{TRUE} to suppress the error.
#' @return the objective function value, or a list of objective function
#' values. They must be numeric and of length one; anything else is an error.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad.stats}}, \code{hydromad.object}
#' @keywords utilities
#' @examples
#'
#' dat <- data.frame(Q = rnorm(10), X = rnorm(10))
#' objFunVal(dat, hmadstat("RMSE"))
#' @export
objFunVal <- function(x, objective, ...) {
  UseMethod("objFunVal")
}


#' @rdname objFunVal
#' @export
objFunVal.default <-
  function(x, objective = hydromad.getOption("objective"),
           ..., nan.ok = FALSE) {
    stopifnot(is.numeric(x) || is.data.frame(x))
    stopifnot(length(colnames(x)) > 0)
    ## these can be referred to in `objective`
    DATA <- x
    X <- x[, "X"]
    delayedAssign("Q", x[, "Q"])
    delayedAssign("U", x[, "U"])
    ## catch the .() function (used for cacheing, see hydromad.stats)
    ## normally it would not get through to here; evaluated by fitBy*() etc.
    ## But this may be needed if objFunVal() is called directly.
    assign(".", function(x) x)
    objFunVal1 <- function(obj, ...) {
      if (inherits(obj, "formula")) {
        val <- with(as.list(environment(obj)), eval(obj[[2]]))
      } else if (is.function(obj)) {
        assign(".", function(x) x, environment(obj))
        val <- obj(Q, X, ..., U = U, DATA = DATA)
      } else {
        stop(
          "'objective' should be a function or formula, not a ",
          toString(class(obj))
        )
      }
      if (is.nan(val)) {
        if (identical(nan.ok, "warn")) {
          warning("objective function returned NaN")
        } else if (!isTRUE(nan.ok)) {
          stop("objective function returned NaN")
        }
      }
      if (!is.numeric(val)) {
        stop("objective should be numeric, not ", toString(class(val)))
      }
      if (length(val) != 1) {
        stop("objective value should be length 1, not ", length(val))
      }
      as.numeric(val)
    }
    if (is.list(objective)) {
      lapply(objective, objFunVal1, ...)
    } else {
      objFunVal1(objective, ...)
    }
  }

## TODO: could this just merge the data and call the default method? slow?
#' @rdname objFunVal
#' @export
objFunVal.hydromad <-
  function(x, objective = hydromad.getOption("objective"),
           ..., all = FALSE, nan.ok = FALSE) {
    model <- x
    ## these can be referred to in `objective`
    X <- fitted(x, all = all)
    if (length(X) == 0) {
      stop("fitted() returned nothing")
    }
    delayedAssign("Q", observed(x, all = all))
    delayedAssign("U", fitted(x, all = all, U = TRUE))
    delayedAssign("DATA", observed(x, all = all, select = TRUE))
    ## catch the .() function (used for cacheing, see hydromad.stats)
    ## normally it would not get through to here; evaluated by fitBy*() etc
    ## But this may be needed if objFunVal() is called directly.
    assign(".", function(x) x)
    isValidModel <- isValidModel(x)
    objFunVal1 <- function(obj, ...) {
      if (!isValidModel) {
        return(NA_real_)
      }
      if (inherits(obj, "formula")) {
        # TODO: avoid converting whole environment to list
        # [names(env) %in% all.vars(obj[[2]])]
        env <- as.list(environment(obj))
        if (any(c("X", "U", "Q", "DATA", "model") %in% names(env))) {
          colliding.vars <- intersect(names(env), c("X", "U", "Q", "DATA", "model"))
          warning(sprintf(
            "Following variables appear to be defined by user: %s\n  objFunVal will use internally calculated values instead.",
            paste(colliding.vars, collapse = ",")
          ))
          env <- env[!names(env) %in% c("X", "U", "Q", "DATA", "model")]
        }
        val <- with(env, eval(obj[[2]]))
      } else if (is.function(obj)) {
        assign(".", function(x) x, environment(obj))
        val <- obj(Q, X, ..., U = U, DATA = DATA, model = model)
      } else {
        stop(
          "'objective' should be a function or formula, not a ",
          toString(class(obj))
        )
      }
      if (is.nan(val)) {
        if (identical(nan.ok, "warn")) {
          warning("objective function returned NaN")
        } else if (!isTRUE(nan.ok)) {
          stop("objective function returned NaN")
        }
      }
      if (!is.numeric(val)) {
        stop("objective should be numeric, not ", toString(class(val)))
      }
      if (length(val) != 1) {
        stop("objective value should be length 1, not ", length(val))
      }
      as.numeric(val)
    }
    if (is.list(objective)) {
      lapply(objective, objFunVal1, ...)
    } else {
      objFunVal1(objective, ...)
    }
  }


objFunVal.tf <- objFunVal.hydromad

#' @rdname objFunVal
#' @export
objFunVal.runlist <- function(x, objective = list(hydromad.getOption("objective"), mean), ...) {
  if (is.list(objective) && length(objective) == 2) {
    switch(hydromad.getOption("parallel")[["objFunVal.runlist"]],
      clusterApply = {
        vals <- parSapply(cl, x, objFunVal, objective[[1]], ...)
      },
      vals <- sapply(x, objFunVal, objective[[1]], ...)
    )
    agg <- objective[[2]](vals)
    stopifnot(length(agg) == 1)
    return(agg)
  } else if (is.function(objective)) {
    return(objective(x, ...))
  }
  stop("Objective is not a list of length 2 or a function")
}

#' @import utils
utils::globalVariables(c("cl"))
