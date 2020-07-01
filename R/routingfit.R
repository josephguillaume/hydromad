## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

doRoutingFit <-
  function(object, inverseFitOnly = FALSE) {
    ## run 'rfit' with appropriate DATA
    routing <- object$routing
    rfit <- object$rfit
    if (is.null(routing)) {
      return(object)
    }
    if (is.null(rfit)) {
      return(object)
    }
    if (!is.list(rfit)) {
      rfit <- list(rfit)
    }
    ## ok, 'rfit' given, fit the routing model.
    ## we should remove any existing routing parameters
    object$parlist <- as.list(coef(object, which = "sma", warn = FALSE, etc = TRUE))
    object$vcov.rfit <- NULL
    doRfit <- function(method, ...) {
      if (missing(method)) {
        stop("missing 'method' in 'rfit' specification")
      }
      if (!is.character(method)) {
        stop("unrecognised 'method' in 'rfit' specification")
      }
      isInverseMethod <- (any(grep("inverse", method)))
      if (isInverseMethod) {
        DATA <- object$data
      } else {
        if (inverseFitOnly) {
          return(NULL)
        }
        DATA <- cbind(
          U = object$U,
          Q = object$data[, "Q"]
        )
      }
      fnName <- paste(routing, method, "fit", sep = ".")
      ## throw an error if we can't find a function of that name
      force(get(fnName, mode = "function"))
      if (isInverseMethod && hydromad.getOption("trace")) {
        message("== fitting routing by ", fnName, " method ==")
      }
      fcall <- quote(RFIT(DATA, ...))
      fcall[[1]] <- as.symbol(fnName)
      ans <- eval(fcall)
      if (isInverseMethod && is.list(ans)) {
        ## these should be set to NULL anyway, but to be sure:
        ans$fitted.values <- NULL
        ans$residuals <- NULL
      }
      ans
    }
    rans <- do.call(doRfit, rfit)
    if (is.null(rans)) {
      return(object)
    }
    if (inherits(rans, "try-error") ||
      inherits(rans, "error")) {
      object$msg <- rans
      return(object)
    }
    ## get fitted parameters; $delay also becomes a parameter
    rcoef <- coef(rans)
    if (!is.null(rans$delay)) {
      rcoef <- c(rcoef, delay = rans$delay)
    }
    object$parlist <- modifyList(object$parlist, as.list(rcoef))
    object$fitted.values <- fitted(rans, all = TRUE)
    ## coerce fitted values back to original data format
    ## (routing fitting functions and tf() coerce to 'ts')
    mostattributes(object$fitted.values) <- attributes(object$data)
    tmp <- try(vcov(rans), silent = TRUE)
    if (inherits(tmp, "try-error")) {
      tmp <- rans$vcov
    }
    object$vcov.rfit <- tmp
    info.rfit <-
      rans[c("warning", "converged", "iteration", "prefilter")]
    object$info.rfit <- info.rfit[!sapply(info.rfit, is.null)]
    ## delete the fit specification, as it has done its work now
    object$used.rfit <- rfit
    object$rfit <- NULL
    ## Add U to data if rfit method wanted it set, accounting for case where object$data had no names
    if (!is.null(rans$U) & is.null(names(object$data))) {
      object$data <- cbind(Q = object$data, U = zoo(rans$U, order.by = index(object$data)))
    } else if (!is.null(rans$U)) {
      object$data <- cbind(object$data, U = zoo(rans$U, order.by = index(object$data)))
    }
    object
  }
