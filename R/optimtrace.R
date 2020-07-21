#' Extract objective function value series from optimisation results
#'
#' Extract objective function value series from optimisation results.
#'
#'
#' @aliases optimtrace optimtrace.default
#' @param object fitted model object. Currently only \code{\link{hydromad}}
#' objects are supported.
#' @param \dots further arguments which may be supported for particular
#' optimisation algorithms. These are: \describe{ \item{list("raw")}{ if
#' \code{TRUE}, then raw function values are returned, which may consist of
#' multiple series in parallel, and may jump up and down. Otherwise (the
#' default), the best value so far at each time step is returned.  }
#' \item{list("objective")}{ a function or formula to calculate from
#' \code{object} to construct the objective function value series.  It is
#' passed to \code{\link{objFunVal}}.  Currently this is only supported for
#' \code{dream} results.  } }
#' @return a \code{zoo} series where the time index represents number of
#' function evaluations.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{fitBySCE}}, \code{\link{fitByDE}},
#' \code{\link{fitByDream}}, \code{\link{fitByOptim}}
#' @keywords optimization
#' @examples
#'
#' ## see examples in fitByOptim, fitBySCE, fitByDE, fitByDream
#' @export
optimtrace <- function(object, ...) {
  UseMethod("optimtrace")
}


#' @export
optimtrace.default <- function(object, ..., model = object) {
  stopifnot(is.list(object))
  stopifnot(length(object$fit.result) > 0)

  getObjSeq(object$fit.result, ..., model = model)
}

getObjSeq <- function(object, ...) {
  UseMethod("getObjSeq")
}

## method for optim() results produced by hydromad
getObjSeq.default <- function(object, ..., raw = FALSE) {
  ppp <- object$objseq
  ans <- zoo(ppp, 1:NROW(ppp))
  if (raw == FALSE) ans <- cummax(na.locf(ans))
  ans
}

## TODO: could re-calculate objective using stored population
getObjSeq.SCEoptim <- function(object, ..., raw = FALSE) {
  funevals <- object$counts
  ppp <- object$POP.FIT.ALL
  if (raw == FALSE) {
    ## best of the chains at each step
    ppp <- do.call("pmin", as.data.frame(ppp))
  }
  ## reverse it to get the maximised objective value
  ppp <- -ppp
  zoo(ppp, seq(0, funevals, length = NROW(ppp) + 1)[-1])
}

## TODO: could re-calculate objective using stored population
getObjSeq.DEoptim <- function(object, ...) {
  funevals <- object$optim$nfeval
  ppp <- object$member$bestvalit
  ## reverse it to get the maximised objective value
  ppp <- -ppp
  zoo(ppp, seq(0, funevals, length = NROW(ppp) + 1)[-1])
}

getObjSeq.dream <-
  function(object, ..., raw = FALSE, objective = NULL, model = NULL) {
    funevals <- object$fun.evals
    ppp <- -object$hist.logp ## negative log likelihood
    if (is.null(objective)) {
      if (raw == FALSE) {
        ## best of the chains at each step
        ppp <- do.call("pmin", c(as.list(as.data.frame(ppp)), na.rm = TRUE))
        ## best result so far at each step
        ppp <- cummin(na.locf(ppp))
      }
    } else {
      ## calculate corresponding objective function values over time.
      stopifnot(!is.null(model))
      objective <- buildCachedObjectiveFun(objective, model)
      ## could do it in a straightforward way by evaluating model at
      ## each time step and for each chain, but would take a long time!
      ## so instead, only re-calculate objective when likelihood improves.
      bestLik <- do.call("pmin", c(as.list(as.data.frame(ppp)), na.rm = TRUE))
      improved <- c(TRUE, diff(cummin(na.locf(bestLik))) != 0)
      ## calculate objective for time steps with improved likelihood
      improv.fval <- sapply(which(improved), function(i) {
        ichain <- which.min(ppp[i, ])
        ## TODO: temporary; don't need this for newly generated results
        ## (there might be multiple minima and some are false, from dream < 0.3-2)
        ichains <- which(ppp[i, ] == min(ppp[i, ]))
        min(sapply(ichains, function(ichain) {
          ipars <- object$Sequences[[ichain]][i, ]
          objFunVal(update(model, newpars = ipars),
            objective = objective
          )
        }))
      })
      fval <- rep(NA_real_, NROW(ppp))
      fval[improved] <- improv.fval
      fval <- na.locf(fval)
      ppp <- fval
    }
    zoo(ppp, seq(0, funevals, length = NROW(ppp) + 1)[-1])
  }
