##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


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
