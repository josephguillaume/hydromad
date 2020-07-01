## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


predict.hydromad <-
  function(object, newdata = NULL,
           which = c("both", "sma", "routing"),
           ..., all = TRUE,
           feasible.set = FALSE,
           glue.quantiles = NULL,
           groups = NULL, FUN = sum,
           return_state = FALSE,
           return_components = FALSE) {
    which <- match.arg(which)
    ## throw an error if any extra arguments
    if (length(list(...)) > 1) {
      stop(
        "unrecognised arguments: ",
        paste(names(list(...)), collapse = ",")
      )
    }
    ## optional aggregation
    doaggr <- identity
    if (!is.null(groups)) {
      doaggr <- function(x) {
        eventapply(x, groups, FUN = FUN)
      }
    }
    if (is.null(newdata)) {
      if (which == "routing") {
        newdata <- object$U
      } else {
        newdata <- object$data
      }
    } else {
      newdata <- as.zooreg(newdata)
    }
    DATA <- newdata
    sma <- object$sma
    routing <- object$routing
    if (which == "routing") sma <- NULL
    if (which == "sma") routing <- NULL
    sma.args <- as.list(coef(object, which = "sma", etc = TRUE))
    r.args <- as.list(coef(object, which = "routing", etc = TRUE))
    if (is.character(sma)) {
      ## construct call to SMA simulation function
      sma.fun <- paste(sma, ".sim", sep = "")
      doSMA <- function(args) {
        ucall <- as.call(c(
          list(
            as.symbol(sma.fun),
            quote(DATA)
          ),
          args
        ))
        if (return_state) {
          ucall$return_state <- TRUE
        }
        ## calculate U
        U <- eval(ucall)
      }
    } else {
      ## default (NULL) SMA action is to return rainfall
      doSMA <- function(args) {
        if (NCOL(DATA) > 1) {
          if ("U" %in% colnames(DATA)) {
            return(DATA[, "U"])
          } else if (("P" %in% colnames(DATA))) {
            return(DATA[, "P"])
          } else {
            stop("No U or P in DATA to be provided to routing")
          }
        } else {
          DATA
        }
      }
    }
    if (is.character(routing)) {
      ## construct call to routing simulation function
      r.fun <- paste(routing, ".sim", sep = "")
      doRouting <- function(U, args) {
        rcall <- as.call(c(
          list(
            as.symbol(r.fun),
            quote(U)
          ),
          args
        ))
        if (return_components) {
          rcall$return_components <- TRUE
        }
        X <- eval(rcall)
      }
    } else {
      ## default (NULL) routing action is to return U (from SMA)
      doRouting <- function(U, args) {
        U
      }
    }
    ## handle full feasible set of parameters -- simple case only
    if (feasible.set) {
      if (is.null(object$feasible.set)) {
        stop("there is no estimate of the feasible set; try defineFeasibleSet()")
      }
      psets <- object$feasible.set
      ## take default arguments from normal fitted coef(); update with psets
      sma.fs.names <- intersect(names(sma.args), colnames(psets))
      r.fs.names <- intersect(names(r.args), colnames(psets))
      ## can use less memory if only estimating max/min bounds
      boundsOnly <- isTRUE(all.equal(glue.quantiles, c(0, 1)))
      sims <- sim.lower <- sim.upper <- NULL
      time <- NULL
      ## TODO: option to show progress bar?
      # result <- lapply(1:NROW(psets), function(i) {
      for (i in 1:NROW(psets)) {
        pset.i <- as.list(psets[i, ])
        ## run SMA
        sma.args.i <- sma.args
        sma.args.i[sma.fs.names] <- pset.i[sma.fs.names]
        U <- doSMA(sma.args.i)
        ## run routing
        r.args.i <- r.args
        r.args.i[r.fs.names] <- pset.i[r.fs.names]
        xsim <- doRouting(U, r.args.i)
        xsim <- doaggr(xsim)
        if (i == 1) {
          time <- index(xsim)
          sim.lower <- sim.upper <- xsim
        }
        if (boundsOnly) {
          sim.lower <- pmin(sim.lower, xsim)
          sim.upper <- pmax(sim.upper, xsim)
        } else {
          sims <- cbind(sims, coredata(xsim))
        }
      }
      if (boundsOnly) {
        ans <- cbind(
          lower = sim.lower,
          upper = sim.upper
        )
      } else if (!is.null(glue.quantiles)) {
        ## weighted quantiles
        glue.quantiles <- sort(glue.quantiles)
        threshold <- object$feasible.threshold
        if (is.infinite(threshold)) {
          threshold <- min(object$feasible.scores, na.rm = TRUE)
        }
        weights <- object$feasible.scores - threshold
        weights <- weights / sum(weights, na.rm = TRUE)
        bounds <- t(apply(sims, 1, safe.wtd.quantile, weights = weights, probs = glue.quantiles, normwt = TRUE))
        colnames(bounds) <- paste("GLUE", glue.quantiles * 100, sep = ".")
        ans <- zoo(bounds, time)
      } else {
        ## glue.quantiles is NULL; return all simulations
        colnames(sims) <- paste("X", 1:NCOL(sims), sep = "")
        ans <- zoo(sims, time)
      }
      return(if (all) ans else stripWarmup(ans, object$warmup))
    }
    ## check that parameters are fully specified
    if (!isFullySpecified(object, which = which)) {
      stop("model parameters are not fully specified")
    }
    ## run SMA
    U <- doSMA(sma.args)
    if (return_state) {
      S <- U
      if (NCOL(S) > 1) {
        stopifnot("U" %in% colnames(S))
        U <- S[, "U"]
      }
    }
    ## run routing
    X <- doRouting(U, r.args)
    if (return_state) {
      if (is.null(routing)) {
        ans <- S
      } else {
        ans <- cbind(S, X)
        if (length(colnames(S)) > 0) {
          colnames(ans)[1:NCOL(S)] <- colnames(S)
        }
      }
    } else {
      ans <- X
    }
    ans <- doaggr(ans)
    if (all) ans else stripWarmup(ans, object$warmup)
  }
