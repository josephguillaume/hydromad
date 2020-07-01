## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

update.hydromad <-
  function(object, ..., newdata = NULL, newpars = NULL,
           sma, routing, rfit, warmup,
           feasible.set, feasible.scores, glue.quantiles,
           and.rescale = TRUE) {
    if (length(newpars) > 0) {
      ## re-call this function with elements of newpars
      if (length(names(newpars)) < length(newpars) ||
        any(names(newpars) == "")) {
        stop("elements of 'newpars' must be named")
      }
      ccall <- match.call()
      ccall[[1]] <- quote(hydromad:::update.hydromad)
      ccall <- as.call(modifyList(
        as.list(ccall),
        as.list(newpars)
      ))
      ccall$newpars <- NULL
      return(eval.parent(ccall))
    }
    ## update timestamp
    object$last.updated <- Sys.time()

    ## update the stored call
    upcall <- match.call()
    nm <- names(upcall)
    if (!is.null(newdata)) {
      object$call$DATA <- upcall["newdata"]
    }
    if (!is.null(nm)) {
      nm <- nm[!(nm %in% c("", "object", "newdata", "and.rescale"))]
      if (length(nm) > 0) {
        extras <- as.list(upcall)[nm]
        ## round off parameter values in call (for display purposes)
        isnum <- sapply(extras, is.numeric)
        extras[isnum] <- lapply(extras[isnum], signif, 6)
        object$call[nm] <- extras
      }
    }
    ## update DATA
    RUNSMA <- is.null(object$U)
    if (!is.null(newdata)) {
      object$data <- as.zooreg(newdata)
      RUNSMA <- TRUE
    }
    ## update named items
    if (!missing(sma)) {
      if (!is.null(object$sma)) {
        ## remove old SMA-specific parameters
        object$parlist <- as.list(coef(object, "routing", etc = TRUE))
        object$feasible.set <- NULL
      }
      class(object) <- "hydromad"
      if (!is.null(sma)) {
        stopifnot(is.character(sma))
        ## take default parameter ranges/values from hydromad.options()
        pardefs <- hydromad.getOption(sma)
        ## or sma.ranges() if available
        if (is.null(pardefs)) {
          pardef.fn <- paste(sma, ".ranges", sep = "")
          if (exists(pardef.fn, mode = "function")) {
            pardefs <- do.call(pardef.fn, list())
          }
        }
        object$parlist <-
          modifyList(
            as.list(pardefs),
            object$parlist
          )
        class(object) <- c(
          paste("hydromad", sma, sep = "."),
          "hydromad"
        )
      }
      object$sma <- sma
      object$sma.fun <- NULL
      object$sma.formals <- NULL
      if (!is.null(sma)) {
        object$sma.fun <- paste(sma, ".sim", sep = "")
        force(get(object$sma.fun, mode = "function")) ## check exists
        object$sma.formals <- formals(object$sma.fun)
      }
      RUNSMA <- TRUE
    }
    if (!missing(routing)) {
      if (!is.null(object$routing)) {
        ## remove old routing-specific parameters
        object$parlist <- as.list(coef(object, "sma", etc = TRUE))
        object$feasible.set <- NULL
      }
      if (!is.null(routing)) {
        stopifnot(is.character(routing))
        ## take default parameter ranges/values from hydromad.options()
        pardefs <- hydromad.getOption(routing)
        ## or routing.ranges() if available
        if (is.null(pardefs)) {
          pardef.fn <- paste(routing, ".ranges", sep = "")
          if (exists(pardef.fn, mode = "function")) {
            pardefs <- do.call(pardef.fn, list())
          }
        }
        object$parlist <-
          modifyList(
            as.list(pardefs),
            object$parlist
          )
      }
      object$routing <- routing
    }
    if (!missing(rfit)) {
      object$rfit <- rfit
    }
    if (!missing(warmup)) {
      stopifnot(is.numeric(warmup))
      object$warmup <- warmup
    }
    ## delete any previous stored results -- we will run again
    object$fitted.values <- NULL
    object$residuals <- NULL
    object$msg <- NULL
    ## allow a fixed routing model to be fit once
    object <- doRoutingFit(object, inverseFitOnly = TRUE)
    ## feasible parameter set (matrix)
    RUNFEASIBLE <- RUNSMA
    if (!missing(feasible.set)) {
      if (!is.null(feasible.set)) {
        stopifnot(is.matrix(feasible.set))
        stopifnot(length(colnames(feasible.set)) > 0)
        RUNFEASIBLE <- TRUE
        if (is.null(object$glue.quantiles)) {
          object$glue.quantiles <- c(0, 1)
        }
      } else {
        object$feasible.fitted <- NULL
      }
      object$feasible.set <- feasible.set
    }
    if (!missing(feasible.scores)) {
      if (!is.null(feasible.scores)) {
        if (length(feasible.scores) == 1) {
          feasible.scores <-
            rep(feasible.scores, NROW(object$feasible.set))
        } else {
          stopifnot(length(feasible.scores) == NROW(object$feasible.set))
        }
        RUNFEASIBLE <- TRUE
      }
      objects$feasible.scores <- feasible.scores
    }
    if (!missing(glue.quantiles)) {
      object$glue.quantiles <- glue.quantiles
      RUNFEASIBLE <- TRUE
    }
    if (RUNFEASIBLE && !is.null(object$feasible.set)) {
      object$feasible.fitted <-
        predict(object,
          feasible.set = TRUE,
          glue.quantiles = object$glue.quantiles
        )
    }
    ## update parameters.
    ## the arguments in `...` may be intended for sma and/or routing
    dots <- list(...)
    ## first, store the SMA parameters to check whether they are changed.
    oldSmaPar <- coef(object, which = "sma", etc = TRUE)
    ## note that any NULL values in dots will delete parlist items
    ## (those parameters will then fall back to sim function defaults)
    object$parlist <- modifyList(object$parlist, dots)
    ## check whether the parameters are fully specified
    if (!isFullySpecified(object, which = "sma")) {
      object$U <- NULL
      return(object)
    }
    ## otherwise... SMA parameters are fully specified.
    ## need to run the SMA if parameters changed.
    newSmaPar <- coef(object, which = "sma", etc = TRUE)
    if (!identical(oldSmaPar, newSmaPar)) {
      RUNSMA <- TRUE
    }
    ## run SMA to generate U
    if (RUNSMA) {
      object$U <- as.zooreg(predict(object,
        which = "sma",
        return_state = isTRUE(object$parlist$return_state)
      ))
    }
    ## handle routing
    routing <- object$routing
    if (!is.null(routing)) {
      ## take starting value of filter from data
      initX <- object$initX
      Xs_0 <- if (is.logical(initX)) 0 else initX
      if (is.null(object$initX)) {
        Q <- observed(object, all = TRUE)
        if (!is.null(Q)) {
          if (any(is.finite(Q[1:10]))) {
            object$initX <- min(Q[1:10], na.rm = TRUE)
          }
        }
      }
      ## TODO: use initX!

      if (!is.null(object$rfit)) {
        ## fit the routing model.
        object <- doRoutingFit(object)
        ## it may have failed:
        if (!isFullySpecified(object)) {
          return(object)
        }
      } else {
        ## just run routing simulation

        ## parameters not fixed, can not run
        if (!isFullySpecified(object)) {
          return(object)
        }

        ## TODO: take starting value of filter from data

        object$fitted.values <-
          as.zooreg(predict(object, which = "routing"))
      }
    }
    ## absorbScale:
    ## move "steady state gain" (volume factor) from unit hydrograph
    ##      -- or model bias if no suitable UH -- into SMA
    ## (if possible; dispatched on SMA class)
    if (!and.rescale) {
      return(object)
    }
    alreadyRescaled <- FALSE
    if (!is.null(routing)) {
      rcoef <- coef(object, which = "routing")
      ssgFn <- paste("ssg", routing, sep = ".")
      if (exists(ssgFn, mode = "function")) {
        gain <- do.call(ssgFn, list(rcoef))
        if (!is.na(gain) && (abs(gain - 1) > 1e-6)) {
          ## UH needs to be normalised (unit volume)
          sobject <- absorbScale(object, gain = gain)
          if (!is.null(sobject)) {
            object <- sobject
          }
          ## Even if we can't absorb scale into SMA:
          ## normalise UH parameters
          normFn <- paste("normalise", routing, sep = ".")
          if (exists(normFn, mode = "function")) {
            normPars <- do.call(normFn, list(rcoef))
            object$parlist <- modifyList(
              object$parlist,
              as.list(normPars)
            )
          } else {
            stop(
              "could not find function ", normFn,
              " (works with ", ssgFn, ")"
            )
          }
          ## note, we might not need to run the model again;
          ## assuming that absorbScale() on SMA has exactly compensated for
          ## normalise() on the routing parameters.
          object <- update(object, and.rescale = FALSE)
          ## NOTE - assuming that after normalise(), ssg.*() will return 1
          ##      - otherwise infinite recursion will result.
          alreadyRescaled <- TRUE
          ## TODO
        }
      }
    }

    ## TODO: i think should leave this out... ?

    if (!alreadyRescaled) {
      ## estimate scale by mass balance
      Q <- observed(object)
      if (!is.null(Q)) {
        X <- fitted(object)
        ok <- complete.cases(Q, X)
        gain <- (sum(Q[ok]) / sum(X[ok]))
        if (!is.finite(gain) || (gain == 0)) {
          return(object)
        }
        sobject <- absorbScale(object, gain = gain)
        if (!is.null(sobject)) {
          object <- sobject
          ## and re-do the routing
          object <- update(object, and.rescale = FALSE)
        }
      }
    }

    return(object)
  }


absorbScale <- function(object, gain, ...) {
  UseMethod("absorbScale")
}

absorbScale.hydromad <- function(object, gain, ...) {
  return(NULL)
}
