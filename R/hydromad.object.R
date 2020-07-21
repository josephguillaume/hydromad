#' Standard methods for Hydromad model objects
#'
#' A \code{hydromad} object represents a model, which may be fully specified
#' (calibrated) or be defined only by parameter ranges.  The model
#' specification and parameter values are stored along with the observed input
#' and output time-series data.
#'
#' Several standard methods are available for \code{hydromad} objects:
#'
#' (note: these are links to the generic functions only)
#'
#' \code{\link{update}}, \code{\link{predict}}, \code{\link{fitted}},
#' \code{\link{observed}}, \code{\link{residuals}}, \code{\link{coef}},
#' \code{\link{vcov}}, etc.
#'
#' The \code{\link[=summary.hydromad]{summary}} and
#' \code{\link[=predict.hydromad]{predict}} methods are documented on different
#' pages.
#'
#' The main plot methods are \code{\link{xyplot.hydromad}} and
#' \code{\link{qqmath.hydromad}}.
#'
#' \code{isValidModel()} returns \code{TRUE} only if the supplied
#' \code{hydromad} object is fully specified and has a calculated output
#' series.
#'
#' To help sample parameters, \code{getFreeParsRanges} returns the list of
#' ranges of parameters that do not have fixed parameter values. In particular,
#' it is used in conjunction with \code{\link{evalPars}} to perform sensitivity
#' analyses.
#'
#' @name hydromad.object
#' @aliases update.hydromad fitted.hydromad observed.hydromad
#' residuals.hydromad coef.hydromad vcov.hydromad
#' isValidModel print.hydromad
#' @param object an object of class \code{hydromad}.
#' @param \dots In the \code{update} method, parameter values or ranges for the
#' SMA and/or routing simulation functions can be given, as with the
#' \code{hydromad()} function.
# @param newdata a \code{\link{ts}}-like object containing a new time series
# dataset (replacing the original \code{DATA} argument given to the
# \code{hydromad} function).
# @param newpars a named list or vector of parameter values; this is
# equivalent to specifying the same values as named arguments (as in
# \dQuote{\dots{}}).
# @param sma,routing,rfit,warmup same arguments as for the
# \code{\link{hydromad}} function. The \code{update} method allows these to be
# changed on an existing model object.
# @param feasible.set,feasible.scores,glue.quantiles the \emph{feasible set}
# of parameter sets can be specified as a matrix, where parameter values are
# given in named columns. The corresponding objective function values for each
# row can be given as \code{feasible.scores}. If \code{glue.quantiles} is
# omitted or NULL, then overall bounds of the ensemble simulation will be
# calculated. Otherwise GLUE-like quantiles can be given as
# \code{glue.quantiles}. See \code{\link{defineFeasibleSet}}.
# @param and.rescale set to \code{FALSE} to suppress any automatic adjustment
# of parameters for mass balance.
# @param which selects either the SMA or routing model, or both models (the
#' default).
#' @param all if \code{TRUE}, return the entire time series for which data
#' exists. Otherwise, the warmup period (specified as an argument to
#' \code{\link{hydromad}} or \code{update}) is stripped off.
#' @param incl.other.vars if \code{TRUE} and model returns a multivariate
#' object (e.g. because of \code{return_components} or \code{return_state}),
#' then return time series for all variables. Otherwise, return only the column
#' named \code{X} or \code{U} (if \code{U=TRUE}). Raises an error if the column
#' is missing.
#' @param feasible.bounds if \code{TRUE}, then ensemble simulation bounds are
#' extracted and returned. This only works if a \emph{feasible set} has been
#' specified using \code{\link{defineFeasibleSet}} or the \code{update} method.
#' Note that the meaning depends on what value of \code{glue.quantiles} was
#' specified to those methods: it might be the overall simulation bounds, or
#' some GLUE-like quantile values. This will be indicated by the returned
#' column names.
#' @param U to return modelled effective rainfall (the output from SMA) rather
#' than streamflow.
#' @param select data series to extract (from the original \code{DATA}
#' argument). Use \code{TRUE} to extract all columns.
# @param warn by default, \code{coef} gives a warning if the model parameters
# are not fully specifed (i.e. some parameters have ranges rather than
# specific values), because it returns a \code{list} rather than a
# \code{vector} in this case. Setting \code{warn = FALSE} skips the warning.
# @param etc by default, \code{coef} returns only the model \emph{parameters},
# which are defined as being numeric and not wrapped in \code{\link{I}()}. If
# \code{etc = TRUE} is given, then all arguments for the simulation
# function(s) will be returned, which may include other data types like
# logicals or time series. In this case the return value is always a
# \code{list}.
#' @param boxcox Placeholder
#' @param start Placeholder
#' @param loglik Placeholder
#' @param x Placeholder
#' @param digits Placeholder
#' @return
#'
#' \code{update} returns a new \code{hydromad} object.
#'
#' \code{fitted}, \code{observed} and \code{residuals} returns time series.
#'
#' \code{coef} returns a named numeric vector, or a named \code{list} if one or
#' more parameters are not fully specified.
#'
#' \code{getFreeParsRanges} returns a named list of parameter ranges, for each
#' parameter that has a range defined. Note that this excludes fixed parameters
#' and parameters defines as a set of discrete values, for which
#' \code{coef(object,warn=FALSE)} should be used instead.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}}, \code{\link{summary.hydromad}},
#' \code{\link{predict.hydromad}}, \code{\link{xyplot.hydromad}},
#' \code{\link{runlist}}, \code{\link{evalPars}}
#' @keywords methods


isFullySpecified <- function(object, ...) {
  !is.list(coef(object, ..., warn = FALSE))
}


#' @rdname hydromad.object
#' @export
fitted.hydromad <-
  function(object, ..., U = FALSE,
           all = FALSE,
           # TODO: change to include.warmup=FALSE - many locations
           feasible.bounds = FALSE,
           incl.other.vars = FALSE) {
    if (is.null(object$routing)) {
      U <- TRUE
    }
    if (!feasible.bounds) {
      # Select either U or X
      tmp <- if (U) object$U else object$fitted.values
      # Return single column from  multivariate objects
      if (is.matrix(tmp) && !incl.other.vars) {
        if (U) {
          if (!"U" %in% names(tmp)) stop("object$U is multivariate and incl.other.vars=F but column U is missing")
          tmp <- tmp[, "U"]
        } else {
          if (!"X" %in% names(tmp)) stop("object$fitted.values is multivariate and incl.other.vars=F but column X is is missing")
          tmp <- tmp[, "X"]
        }
      }
    } else if (feasible.bounds) {
      if (is.null(object$feasible.fitted)) {
        stop("there is no estimate of the feasible bounds; try defineFeasibleSet()")
      }
      tmp <- object$feasible.fitted
    }
    if (length(tmp) == 0) {
      return(tmp)
    }
    if (all) tmp else stripWarmup(tmp, object$warmup)
  }


#' @rdname hydromad.object
#' @export
residuals.hydromad <-
  function(object, ..., all = FALSE, boxcox = FALSE, start = NULL) {
    fit <- fitted(object, all = TRUE)
    if (length(fit) == 0) {
      return(fit)
    }
    obs <- object$data[, "Q"]
    if (!identical(boxcox, FALSE)) {
      coreQ <- coredata(na.omit(obs))
      if (is.null(start)) {
        start <-
          quantile(coreQ[coreQ > 0], 0.1, names = FALSE)
      }
      if (isTRUE(boxcox)) {
        lambda <- coef(powerTransform(coreQ + start))
      } else if (is.numeric(boxcox)) {
        lambda <- boxcox
      } else {
        stop("'boxcox' should be logical or numeric")
      }
      coredata(obs) <- bcPower(coredata(obs) + start, lambda)
      coredata(fit) <- bcPower(coredata(fit) + start, lambda)
    }
    tmp <- (obs - fit)
    if (all) tmp else stripWarmup(tmp, object$warmup)
  }


#' @rdname hydromad.object
#' @export
observed.hydromad <- function(object, ..., select = "Q", all = FALSE) {
  ## observed.default will work (for Q), but this may be slightly faster
  if (is.character(select)) {
    if (!all(select %in% colnames(object$data))) {
      return(NULL)
    }
  }
  tmp <- object$data[, select]
  if (all) tmp else stripWarmup(tmp, object$warmup)
}


#' @rdname hydromad.object
#' @export
vcov.hydromad <- function(object, ...) {
  cov.mat <- object$cov.mat ## may be NULL
  rcov <- object$vcov.rfit
  ans <- cov.mat
  if (!is.null(rcov)) {
    if (is.null(ans)) {
      ans <- rcov
    } else {
      ## merge the two matrices
      tmp.right <- rbind(
        matrix(ncol = ncol(rcov), nrow = nrow(cov.mat)),
        rcov
      )
      ans <- rbind(cov.mat, matrix(ncol = ncol(cov.mat), nrow = nrow(rcov)))
      ans <- cbind(ans, tmp.right)
      rownames(ans) <- colnames(ans) <-
        c(rownames(cov.mat), rownames(rcov))
    }
  }
  ans
}


#' @rdname hydromad.object
#' @export
logLik.hydromad <-
  function(object, loglik = hydromad.getOption("loglik"), ...) {
    val <- objFunVal(object, objective = loglik, ...)
    ## TODO: for a fitted model we do not know how many parameters were fitted
    ## guess:
    attr(val, "df") <- length(coef(object))
    attr(val, "nobs") <- attr(val, "nall") <- length(fitted(object, all = TRUE))
    class(val) <- "logLik"
    val
  }


#' @rdname hydromad.object
#' @export
deviance.hydromad <- stats:::deviance.lm


#' @rdname hydromad.object
#' @export
print.hydromad <-
  function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\n",
      "Hydromad model with ", toString(deparse(x$sma)), " SMA",
      " and ", toString(deparse(x$routing)), " routing:", "\n",
      sep = ""
    )
    rx <- x$data
    cat("Start = ", index2char(index(rx)[1], frequency(rx)),
      ", End = ", index2char(index(rx)[NROW(rx)], frequency(rx)),
      "\n",
      sep = ""
    )
    cat("\n")
    for (which in c("sma", "routing")) {
      if (!is.null(x[[which]])) {
        if (which == "sma") {
          cat("SMA Parameters:\n")
        } else {
          cat("Routing Parameters:\n")
        }
        if (isFullySpecified(x, which = which)) {
          ## all unique parameter values
          coefx <- coef(x, which = which)
          print(coefx, digits = digits, quote = FALSE, print.gap = 2)
          if ((which == "routing") &&
            isTRUE(x$routing %in% c("armax", "expuh"))) {
            tmp <- describeTF(coefx)
            if (!is.null(tmp) && !is.na(tmp)) {
              cat("TF Structure:", tmp, "\n")
            }
          }
        } else {
          ## one or more parameters specified as ranges only
          parlist <- coef(x, which = which, warn = FALSE)
          print(data.frame(
            lower = sapply(parlist, min),
            upper = sapply(parlist, max),
            ` ` = sapply(parlist, function(p) {
              ## mark fixed parameters
              if (isTRUE(diff(range(p)) == 0)) "(==)" else ""
            }), check.names = FALSE
          ),
          digits = digits
          )
        }
        # cat("\n")
      }
    }
    if (!is.null(x$feasible.set)) {
      cat("Feasible parameter set:\n")
      print(apply(x$feasible.set, 2, function(xi) {
        signif(c(lower = min(xi), upper = max(xi)), digits = digits)
      }))
    }
    if (!is.null(x$rfit)) {
      cat(
        "Routing fit spec.:",
        toString(deparse(x$rfit, control = c(), width = 500),
          width = getOption("width")
        ), "\n"
      )
    }
    if (!is.null(x$fit.call)) {
      cat("\nFit: ($fit.result)\n")
      print(x$fit.call)
      cat(
        "    ", x$funevals, "function evaluations in",
        x$timing[3], "seconds", "\n"
      )
    }
    if (length(x$info.rfit) > 0) {
      cat(
        "\nRouting fit info: ",
        toString(deparse(x$info.rfit, control = c(), width = 500),
          width = getOption("width")
        ), "\n"
      )
    }
    if (!is.null(x$msg)) {
      cat("\nMessage:", toString(x$msg), "\n")
    }
    invisible(x)
  }


#' @rdname hydromad.object
#' @export
str.hydromad.runlist <-
  function(object, ...) {
    cat("\nList of Hydromad model runs:\n")
    str(lapply(object, function(obj) {
      if (!is.null(obj$msg)) {
        list(call = obj$call, message = obj$msg)
      } else {
        obj$call
      }
    }))
    invisible()
  }


#' @rdname hydromad.object
#' @export
isValidModel <- function(object, ...) {
  UseMethod("isValidModel")
}


#' @rdname hydromad.object
#' @export
isValidModel.default <- function(object, ...) {
  return(FALSE)
}


#' @rdname hydromad.object
#' @export
isValidModel.hydromad <- function(object, ...) {
  is.numeric(fitted(object, all = TRUE))
}
