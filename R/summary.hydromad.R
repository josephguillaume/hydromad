## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

#' Assess and summarise performance of Hydromad models
#'
#' Assess and summarise performance of Hydromad models.
#'
#' Definitions of statistics are looked up by name in
#' \code{\link{hydromad.stats}()}, and you can also use that function to add
#' new statistics.
#'
#' @importFrom zoo time<- coredata<-
#'
#' @aliases summary.hydromad print.summary.hydromad
#' print.summaryWithBreaks.hydromad summary.hydromad.runlist
#' @param object an object of class \code{hydromad}.
#' @param breaks if specified, break up the time series and calculate
#' statistics for each chunk.  Can be a vector of cut points \emph{or} number
#' of intervals to cut into \emph{or} an interval specification, one of "sec",
#' "min", "hour", "day", "DSTday", "week", "month" or "year", optionally
#' preceded by an integer and a space, or followed by "s".  See
#' \code{\link{cut.Date}}.
#' @param stats which performance statistics to compute. The definitions are
#' looked up by name in \code{\link{hydromad.stats}()}. See Details.
#' @param with.hydrostats to also include \code{timesteps}, \code{missing}
#' (number of timesteps with missing data), \code{mean.P}, \code{mean.Q},
#' \code{runoff} (mean rainfall and streamflow, and their ratio).
#' @param na.action an optional function to apply to the data to treat missing
#' values before calculating statistics.
#' @param \dots ignored.
#' @return \code{summary} returns a list with named entries for each of the
#' chosen stats (\code{stats}). When \code{breaks} is given it retuns a
#' \code{zoo} object of class \code{"summaryWithBreaks"}.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad.stats}}, \code{\link{objFunVal}} for a way to
#' calculate statistic values directly; \code{hydromad.object}
#' @keywords methods
#' @examples
#'
#'
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData,
#'   sma = "scalar",
#'   routing = "expuh", tau_s = 10
#' )
#' mod0
#'
#' summary(mod0)
#'
#' summary(mod0, breaks = "months")
#' # summary(mod0, breaks = 3)
#'
#' allstats <- names(hydromad.stats())
#' ## Ignore r.sq.vartd because it requires event to be specified
#' allstats <- setdiff(allstats, "r.sq.vartd")
#' summary(mod0, stats = allstats)
#'
#' objFunVal(mod0, hmadstat("r.sq.log"))
#'
#' hydromad_stats <- hydromad.stats()
#' # Ignore r.sq.vartd because it requires event to be specified
#' hydromad_stats$r.sq.vartd <- NULL
#' str(objFunVal(mod0, hydromad_stats))
#' @export
summary.hydromad <-
  function(object, breaks = NULL,
           stats = hydromad.getOption("summary.stats"),
           with.hydrostats = TRUE,
           na.action = na.exclude,
           ...) {
    if (!isFullySpecified(object)) {
      stop("model parameters not fully specified")
    }
    if (!isValidModel(object)) {
      stop("model is not valid")
    }

    hydrostats <- function(chunk) {
      ok <- complete.cases(chunk)
      meanP <- mean(chunk[ok, "P"])
      meanQ <- mean(chunk[ok, "Q"])
      list(
        timesteps = NROW(chunk),
        missing = sum(!ok),
        mean.P = meanP,
        mean.Q = meanQ,
        runoff = meanQ / meanP
      )
    }

    ## pull out definitions of statistics
    stats2 <- setdiff(stats, c("ARPE", "YIC"))
    stats.def <- hydromad.stats(stats2)
    bad <- sapply(stats.def, is.null)
    if (any(bad)) {
      stop(
        "no definition found in hydromad.stats() for: ",
        toString(stats[bad])
      )
    }

    DATA <- cbind(
      P = observed(object, select = "P", all = TRUE),
      Q = observed(object, all = TRUE),
      X = fitted(object, all = TRUE),
      U = fitted(object, U = TRUE, all = TRUE)
    )

    if (is.null(breaks)) {
      ## remove warmup
      DATA <- stripWarmup(DATA, object$warmup)
    } else {
      group <- cut(time(DATA), breaks = breaks)
      group <- factor(group)
      ## remove warmup
      DATA <- stripWarmup(DATA, object$warmup)
      if (object$warmup > 0) {
        group <- group[-seq_len(object$warmup)]
      }
      ## remove any groups with only one data point
      bad <- table(group) <= 1
      levels(group)[bad] <- NA

      ans <-
        eventapply(DATA, group,
          FUN = function(x) {
            unlist(c(
              if (with.hydrostats) hydrostats(x),
              objFunVal(x,
                objective = stats.def,
                nan.ok = "warn", model = object
              )
            ))
          },
          by.column = FALSE
        )
      ## copy the last entry with the final date, to mark the end of last period
      lastbit <- tail(ans, 1)
      time(lastbit) <- end(DATA)
      ans <- rbind(ans, lastbit)
      return(ans)
    }

    ans <- list(call = object$call)

    ## ARPE and YIC
    if (any(c("ARPE", "YIC") %in% stats)) {
      arpe <- NA_real_
      if (!is.null(vcov(object))) {
        coef.var <- diag(vcov(object))
        cc <- coef(object, "routing")
        cc2 <- tfParsConvert(cc, "a,b")
        cc[names(cc2)] <- cc2
        nms <- intersect(names(cc), names(coef.var))
        arpe <- mean(coef.var[nms] / (cc[nms]^2))
      }
      if ("ARPE" %in% stats) {
        ans$ARPE <- arpe
      }
      if ("YIC" %in% stats) {
        var.ratio <- (var(residuals(object), na.rm = TRUE) /
          var(observed(object), na.rm = TRUE))
        ans$YIC <- log(var.ratio) + log(arpe)
      }
    }

    if (with.hydrostats) {
      ans <- c(ans, hydrostats(DATA))
    }

    ## call objFunVal for the rest
    ans <- c(ans, objFunVal(object, objective = stats.def, nan.ok = "warn"))

    class(ans) <- "summary.hydromad"
    ans
  }


#' @export
print.summary.hydromad <-
  function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    if (!is.null(x$timesteps) && !is.null(x$missing)) {
      cat("Time steps: ", x$timesteps, " (", x$missing, " missing).\n", sep = "")
    }
    if (!is.null(x$mean.P) && !is.null(x$mean.Q)) {
      cat("Runoff ratio (Q/P): (",
        format(x$mean.Q, digits = digits), " / ",
        format(x$mean.P, digits = digits), ") = ",
        format(x$mean.Q / x$mean.P, digits = digits), "\n",
        sep = ""
      )
    }
    if (!is.null(x$YIC)) {
      cat("YIC:", format(x$YIC, digits = digits), "\n")
    }
    if (!is.null(x$ARPE)) {
      cat("ARPE (%):", format(x$ARPE * 100, digits = digits), "\n")
    }
    ## remove these already-shown ones
    nms <- setdiff(
      names(x),
      c(
        "call", "timesteps", "missing",
        "mean.Q", "mean.P", "runoff",
        "YIC", "ARPE"
      )
    )
    for (nm in nms) {
      xi <- x[[nm]]
      if (is.numeric(xi) && length(xi) == 1) {
        nm <- gsub("\\.", " ", nm)
        cat(nm, ": ", format(xi, digits = digits), "\n", sep = "")
      }
    }
    cat("\n", "For definitions see ?hydromad.stats\n", "\n", sep = "")
    invisible(x)
  }


#' @export
print.summaryWithBreaks.hydromad <-
  function(x, digits = max(3, getOption("digits") - 3), ...) {
    ## just simplify the printed output by rounding
    NextMethod("print")
  }


#' @export
summary.hydromad.runlist <-
  function(object, with.hydrostats = FALSE, ...) {
    summary.runlist(object, ...,
      with.hydrostats = with.hydrostats
    )
  }
