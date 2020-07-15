## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##



#' Specify rainfall - runoff (hydrology) models.
#' 
#' The \code{hydromad} function can be used to specify models with their model
#' equations, data, parameters and settings. It allows a general two-component
#' structure, where the Soil Moisture Accounting (\code{sma}) component and the
#' Routing (\code{routing}) component can be arbitrary functions. A method can
#' be specified for fitting the dependent routing component.
#' 
#' The \code{hydromad()} function allows models to be specified with the given
#' component models and parameter specifications. The resulting object can
#' later be modified using the \code{update.hydromad} method
#' using the same syntax.
#' 
#' Methods for working with the model objects are listed under
#' \code{hydromad.object}.
#' 
#' For a tutorial, type \code{vignette("tutorial", package = "hydromad")}.
#' 
#' For an overview of the package, see the paper
#' \code{vignette("hydromad_paper")}.
#' 
#' For a list of the package functions with their help pages, see the website
#' \url{http://hydromad.catchment.org/}.
#' 
#' @param DATA a \code{\link{ts}}-like object with named columns: \describe{
#' \item{list("P")}{ time series of areal rainfall depths, usually in mm.  }
#' \item{list("E")}{ (optional) time series of potential evapo-transpiration,
#' or more typically, temperature as an indicator of this. Required for some
#' models but not others.  } \item{list("Q")}{ (optional) time series of
#' discharge (streamflow) at the catchment outlet. Required for calibration but
#' not simulation.  It should usually be in units of mm (averaged over the
#' catchment area). Use \code{\link{convertFlow}} to convert it.  }
#' \item{etc.}{ other data columns may also be included, and will be accessible
#' via the \code{observed()} method.  } }
#' @param \dots values or ranges for named parameters. Any parameters not given
#' here will be taken from defaults given in \code{hydromad.options(sma)}
#' and/or \code{hydromad.options(routing)}. In addition, other arbitrary
#' arguments may be given here that will be passed on to the simulation
#' function(s) and not treated as parameters. To specify a numeric object that
#' is not a parameter (such as a time series object), wrap it in
#' \code{\link{I}()}.
#' @param sma name of the Soil Moisture Accounting (SMA) component. May be
#' \code{NULL}, in which case the input rainfall will be passed directly to
#' \code{routing}. If \code{sma} is specified, a corresponding simulation
#' function \var{sma}\code{.sim} must exist.
#' @param routing name of the routing component (i.e. the component which takes
#' in effective rainfall from \code{sma} and converts it to streamflow).  May
#' be \code{NULL}, in which case the model output is taken as the output from
#' \code{sma} directly.
#' @param rfit optional specification for fitting the routing component. If a
#' character string is given, then a corresponding function
#' \var{routing}\code{.}\var{rfit}\code{.fit} must exist.
#' @param warmup warmup period in number of time steps.
#' @return the result from \code{hydromad()} is a
#' \code{hydromad object}.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{hydromad.object}
#' @references F.T. Andrews, B.F.W. Croke and A.J. Jakeman (2011). An open
#' software environment for hydrological model assessment and development.
#' \emph{Environmental Modelling and Software} 26 (2011), pp. 1171-1185.
#' \url{http://dx.doi.org/10.1016/j.envsoft.2011.04.006}
#' @keywords models
#' @examples
#' data(Cotter)
#' x <- Cotter[1:1000]
#' 
#' ## IHACRES CWI model with exponential unit hydrograph
#' ## an unfitted model, with ranges of possible parameter values
#' modx <- hydromad(x, sma = "cwi", routing = "expuh",
#'                  tau_s = c(2,100), v_s = c(0,1))
#' modx
#' ## now try to fit it
#' fitx <- fitByOptim(modx)
#' fitx
#' summary(fitx)
#' xyplot(fitx, with.P = TRUE, type = c("l", "g"))
#' 
#' data(Canning)
#' x <- window(Canning, start = "1980-01-01", end = "1982-01-01")
#' xyplot(x)
#' ## IHACRES CWI model with extra parameter l
#' ## Fixed UH (fit once) by inverse method
#' ## an unfitted model, with ranges of possible parameter values
#' mod0 <- hydromad(x, sma = "cwi", l = c(0, 100), 
#'              routing = "armax", rfit = list("inverse", order = c(1,1)))
#' mod0
#' ## now try to fit the free parameters
#' fit1 <- fitByOptim(mod0)
#' fit1
#' summary(fit1)
#' xyplot(fit1)
#' 
#' @export hydromad
hydromad <-
  function(DATA = zoo(),
           ...,
           sma = hydromad.getOption("sma"),
           routing = hydromad.getOption("routing"),
           rfit = NULL,
           warmup = hydromad.getOption("warmup")) {
    ## create the model object
    obj <- list(call = match.call())
    class(obj) <- "hydromad"
    ## dots `...` may contain arguments for sma and/or routing.
    ## update() takes default parameter ranges/values from hydromad.options().
    obj$parlist <- list()
    obj <- update(obj, ...,
      newdata = DATA, sma = sma,
      routing = routing, rfit = rfit,
      warmup = warmup
    )
    obj$call <- match.call() ## reset call after update()
    return(obj)
  }

isFullySpecified <- function(object, ...) {
  !is.list(coef(object, ..., warn = FALSE))
}

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

deviance.hydromad <- stats:::deviance.lm

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

isValidModel <- function(object, ...) {
  UseMethod("isValidModel")
}

isValidModel.default <- function(object, ...) {
  return(FALSE)
}

isValidModel.hydromad <- function(object, ...) {
  is.numeric(fitted(object, all = TRUE))
}
