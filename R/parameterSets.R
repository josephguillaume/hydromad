##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##




#' Generate parameter sets
#'
#' Generate parameter sets from given ranges, with chosen sampling scheme.
#'
#' Method \code{"latin.hypercube"} generates a regular sequence for each free
#' parameter range (using \code{\link{quantile}} and \code{\link{ppoints}}),
#' and a repeated sequence for each fixed parameter set. Each of these
#' sequences is then randomly permuted. For the special case of \code{samples =
#' 1}, the mean of each range is returned.
#'
#' Method \code{"random"} generates uniform random values in each free
#' parameter range, and random samples from each fixed parameter set.
#'
#' Method \code{"all.combinations"} generates a regular sequence for each free
#' parameter range, and keeps each fixed parameter set as given. All
#' combinations of these values are then calculated (using \code{expand.grid}).
#' The length of the free parameter sequences is chosen such that the total
#' number of results does not exceed \code{samples}. However, if fixed
#' parameter sets are given, the number of combinations of these may exceed
#' \code{samples}, and then any free parameters will be fixed at their mean
#' values.
#'
#' Replicable results can be obtained by using \code{\link{set.seed}}.
#'
#' @param par.ranges A named list of (numeric) parameter values. Each element
#' can be: \itemize{ \item a single number, for fixed parameters; \item a
#' length-2 vector, representing a range of values; \item a vector of length >
#' 2, representing a fixed set of values. Also, a length 2 vector can be
#' defined as a fixed set of values by wrapping it in \code{I()}. Such fixed
#' sets of values are resampled when \code{method = "latin.hypercube"} or
#' \code{"random"}.  }
#' @param samples Number of samples to generate. In the
#' \code{"all.combinations"} method, the result may not have exactly this
#' length.
#' @param method the sampling scheme; see Details.
#' @return the result is a \code{data.frame}, with named columns for each
#' parameter in \code{par.ranges}. Each row represents one parameter set.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{sample}}, \code{\link{quantile}},
#' \code{\link{evalPars}} to then evaluate the parameter sets with a model in
#' hydromad
#' @keywords utilities
#' @examples
#'
#' pars <- list(a = 1, b = 0:1, c = c(1, 10), d = c(-2, -4, -8))
#'
#' set.seed(10)
#' hydromad::parameterSets(pars, 10, method = "random")
#'
#' hydromad::parameterSets(pars, 10, method = "latin")
#'
#' hydromad::parameterSets(pars, 10, method = "all.combinations")
#'
#' hydromad::parameterSets(pars, 20, method = "all.combinations")
#' @export
parameterSets <-
  function(par.ranges, samples,
           method = c("latin.hypercube", "random", "all.combinations")) {
    method <- match.arg(method)
    stopifnot(is.list(par.ranges))
    stopifnot(all(sapply(par.ranges, length) > 0))
    stopifnot(is.numeric(samples))
    ## find parameters that have a null (zero) range, i.e. fixed
    fixed <- sapply(par.ranges, function(x) diff(range(x)) == 0)
    ## find parameters with two-element range, interpreted as a free range
    free <- sapply(par.ranges, function(x) {
      !inherits(x, "AsIs") && length(x) == 2 && (diff(range(x)) > 0)
    })
    ## the rest have length > 2 or AsIs, interpreted as a specified set of values
    spec <- !free & !fixed

    if (method == "all.combinations") {
      ## work out number of samples for free params in order to keep under total 'samples'
      freesamp <- NA
      if (any(free)) {
        specsamp <- prod(unlist(lapply(par.ranges[spec], length)))
        freesamp <- 1
        repeat {
          if (specsamp * ((freesamp + 1)^sum(free)) > samples) {
            break
          }
          freesamp <- freesamp + 1
        }
      }
      par.seqs <- list()
      for (p in names(par.ranges)) {
        vv <- par.ranges[[p]]
        if (fixed[[p]]) {
          ## if parameter is fixed, leave as one value
          par.seqs[[p]] <- vv[1]
        } else if (free[[p]]) {
          if (samples == 1) {
            ## special case of one sample
            par.seqs[[p]] <- mean(vv)
          } else {
            par.seqs[[p]] <-
              # zapsmall(quantile(vv, ppoints(freesamp), names = FALSE))
              zapsmall(seq(min(vv), max(vv), length = freesamp))
          }
        } else {
          par.seqs[[p]] <- vv
        }
      }
      psets <- expand.grid(par.seqs, KEEP.OUT.ATTRS = FALSE)
    }

    if (method == "latin.hypercube") {
      par.seqs <- list()
      for (p in names(par.ranges)) {
        vv <- par.ranges[[p]]
        if (samples == 1) {
          ## special case of one sample
          par.seqs[[p]] <- mean(vv)
        } else if (free[[p]]) {
          par.seqs[[p]] <-
            # zapsmall(quantile(vv, ppoints(samples), names = FALSE))
            zapsmall(seq(min(vv), max(vv), length = samples))
        } else {
          par.seqs[[p]] <- rep(vv, length = samples)
        }
      }
      psets <- data.frame(par.seqs)
      if (samples > 1) {
        psets <- data.frame(lapply(psets, sample))
      }
    }

    if (method == "random") {
      par.seqs <- list()
      for (p in names(par.ranges)) {
        vv <- as.numeric(par.ranges[[p]])
        if (free[[p]]) {
          par.seqs[[p]] <-
            zapsmall(runif(samples, min = min(vv), max = max(vv)))
        } else if (length(vv) == 1) {
          par.seqs[[p]] <- vv
        } else {
          par.seqs[[p]] <-
            sample(vv, samples, replace = TRUE)
        }
      }
      psets <- data.frame(par.seqs)
    }

    return(psets)
  }
