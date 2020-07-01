
## Translated from MATLAB to R
## and substantially revised by Felix Andrews <felix@nfrac.org>
## 2009-08-18

## Changed sampling scheme of parents from each complex;
## convergence criteria; memory efficiency; initial sampling; etc.


# Copyright (C) 2006 Brecht Donckels, BIOMATH, brecht.donckels@ugent.be
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
# USA.


## Originally based on code by Qingyun Duan, 16 May 2005
## http://www.mathworks.com.au/matlabcentral/fileexchange/7671


## NOTE: keep this in sync with the help page!
sceDefaults <- function() {
  list(
    ncomplex = 5, ## number of complexes
    cce.iter = NA, ## number of iteration in inner loop (CCE algorithm)
    fnscale = 1, ## function scaling factor (set to -1 for maximisation)
    elitism = 1, ## controls amount of weighting in sampling towards the better parameter sets
    initsample = "latin", ## sampling scheme for initial values -- "latin" or "random"
    reltol = 1e-5, ## convergence threshold: relative improvement factor required in an SCE iteration
    tolsteps = 7, ## number of iterations within reltol to confirm convergence
    maxit = 10000, ## maximum number of iterations
    maxeval = Inf, ## maximum number of function evaluations
    maxtime = Inf, ## maximum duration of optimization in seconds
    returnpop = FALSE, ## whether to return populations from all iterations
    trace = 0, ## level of user feedback
    REPORT = 1
  )
} ## number of iterations between reports when trace >= 1


SCEoptim <- function(FUN, par, ...,
                     lower = -Inf, upper = Inf,
                     control = list()) {
  FUN <- match.fun(FUN)
  stopifnot(is.numeric(par))
  stopifnot(length(par) > 0)
  stopifnot(is.numeric(lower))
  stopifnot(is.numeric(upper))
  ## allow `lower` or `upper` to apply to all parameters
  if (length(lower) == 1) {
    lower <- rep(lower, length = length(par))
  }
  if (length(upper) == 1) {
    upper <- rep(upper, length = length(par))
  }
  stopifnot(length(lower) == length(par))
  stopifnot(length(upper) == length(par))

  ## determine number of variables to be optimized
  NDIM <- length(par)

  ## update default options with supplied options
  stopifnot(is.list(control))
  control <- modifyList(sceDefaults(), control)
  isValid <- names(control) %in% names(sceDefaults())
  if (any(!isValid)) {
    stop(
      "unrecognised options: ",
      toString(names(control)[!isValid])
    )
  }

  returnpop <- control$returnpop
  trace <- control$trace
  nCOMPLEXES <- control$ncomplex
  CCEITER <- control$cce.iter
  MAXIT <- control$maxit
  MAXEVAL <- control$maxeval

  ## recommended number of CCE steps in Duan et al 1994:
  if (is.na(CCEITER)) {
    CCEITER <- 2 * NDIM + 1
  }

  if (is.finite(MAXEVAL)) {
    ## upper bound on number of iterations to reach MAXEVAL
    MAXIT <- min(MAXIT, ceiling(MAXEVAL / (nCOMPLEXES * CCEITER)))
  }

  ## define number of points in each complex
  nPOINTS_COMPLEX <- 2 * NDIM + 1

  ## define number of points in each simplex
  nPOINTS_SIMPLEX <- NDIM + 1

  ## define total number of points
  nPOINTS <- nCOMPLEXES * nPOINTS_COMPLEX

  ## initialize counters
  funevals <- 0


  costFunction <- function(FUN, par, ...) {
    ## check lower and upper bounds
    i <- which(par < lower)
    if (any(i)) {
      i <- i[1]
      return(1e12 + (lower[i] - par[i]) * 1e6)
    }
    i <- which(par > upper)
    if (any(i)) {
      i <- i[1]
      return(1e12 + (par[i] - upper[i]) * 1e6)
    }
    funevals <<- funevals + 1
    result <- FUN(par, ...) * control$fnscale
    if (is.na(result)) {
      result <- 1e12
    }
    result
  }

  simplexStep <- function(P, FAC) {
    ## Extrapolates by a factor FAC through the face of the simplex across from
    ## the highest (i.e. worst) point.
    worst <- nPOINTS_SIMPLEX
    centr <- apply(P[-worst, , drop = FALSE], 2, mean)
    newpar <- centr * (1 - FAC) + P[worst, ] * FAC
    newpar
  }


  ## initialize population matrix
  POPULATION <- matrix(as.numeric(NA), nrow = nPOINTS, ncol = NDIM)
  if (!is.null(names(par))) {
    colnames(POPULATION) <- names(par)
  }
  POP.FITNESS <- numeric(length = nPOINTS)
  POPULATION[1, ] <- par

  ## generate initial parameter values by random uniform sampling
  finitelower <- ifelse(is.infinite(lower), -(abs(par) + 2) * 5, lower)
  finiteupper <- ifelse(is.infinite(upper), +(abs(par) + 2) * 5, upper)
  if (control$initsample == "latin") {
    for (i in 1:NDIM) {
      tmp <- seq(finitelower[i], finiteupper[i], length = nPOINTS - 1)
      tmp <- jitter(tmp, factor = 2)
      tmp <- pmax(finitelower[i], pmin(finiteupper[i], tmp))
      POPULATION[-1, i] <- sample(tmp)
    }
  } else {
    for (i in 1:NDIM) {
      POPULATION[-1, i] <- runif(nPOINTS - 1, finitelower[i], finiteupper[i])
    }
  }

  ## only store all iterations if requested -- could be big!
  if (!is.finite(MAXIT)) {
    MAXIT <- 10000
    warning("setting maximum iterations to 10000")
  }
  if (returnpop) {
    POP.ALL <- array(as.numeric(NA), dim = c(nPOINTS, NDIM, MAXIT))
    if (!is.null(names(par))) {
      dimnames(POP.ALL)[[2]] <- names(par)
    }
  }
  POP.FIT.ALL <- matrix(as.numeric(NA), ncol = nPOINTS, nrow = MAXIT)
  BESTMEM.ALL <- matrix(as.numeric(NA), ncol = NDIM, nrow = MAXIT)
  if (!is.null(names(par))) {
    colnames(BESTMEM.ALL) <- names(par)
  }

  ## the output object
  obj <- list()
  class(obj) <- c("SCEoptim", class(obj))
  obj$call <- match.call()
  obj$control <- control

  EXITFLAG <- NA
  EXITMSG <- NULL

  ## initialize timer
  tic <- as.numeric(Sys.time())
  toc <- 0

  ## calculate cost for each point in initial population
  for (i in 1:nPOINTS) {
    POP.FITNESS[i] <- costFunction(FUN, POPULATION[i, ], ...)
  }

  ## sort the population in order of increasing function values
  idx <- order(POP.FITNESS)
  POP.FITNESS <- POP.FITNESS[idx]
  POPULATION <- POPULATION[idx, , drop = FALSE]

  ## store one previous iteration only
  POP.PREV <- POPULATION
  POP.FIT.PREV <- POP.FITNESS

  if (returnpop) {
    POP.ALL[, , 1] <- POPULATION
  }
  POP.FIT.ALL[1, ] <- POP.FITNESS
  BESTMEM.ALL[1, ] <- POPULATION[1, ]

  ## store best solution from last two iterations
  prevBestVals <- rep(Inf, control$tolsteps)
  prevBestVals[1] <- POP.FITNESS[1]

  ## for each iteration...
  i <- 0
  while (i < MAXIT) {
    i <- i + 1

    ## The population matrix POPULATION will now be rearranged into complexes.

    ## For each complex ...
    for (j in 1:nCOMPLEXES) {

      ## construct j-th complex from POPULATION

      k1 <- 1:nPOINTS_COMPLEX
      k2 <- (k1 - 1) * nCOMPLEXES + j

      COMPLEX <- POP.PREV[k2, , drop = FALSE]
      COMPLEX_FITNESS <- POP.FIT.PREV[k2]

      ## Each complex evolves a number of steps according to the competitive
      ## complex evolution (CCE) algorithm as described in Duan et al. (1992).
      ## Therefore, a number of 'parents' are selected from each complex which
      ## form a simplex. The selection of the parents is done so that the better
      ## points in the complex have a higher probability to be selected as a
      ## parent. The paper of Duan et al. (1992) describes how a trapezoidal
      ## probability distribution can be used for this purpose.

      for (k in 1:CCEITER) {

        ## select simplex by sampling the complex

        ## sample points with "trapezoidal" i.e. linear probability
        weights <- rev(ppoints(nPOINTS_COMPLEX))
        ## 'elitism' parameter can give more weight to the better results:
        weights <- weights^control$elitism
        LOCATION <- sample(seq(1, nPOINTS_COMPLEX),
          size = nPOINTS_SIMPLEX,
          prob = weights
        )

        LOCATION <- sort(LOCATION)

        ## construct the simplex
        SIMPLEX <- COMPLEX[LOCATION, , drop = FALSE]
        SIMPLEX_FITNESS <- COMPLEX_FITNESS[LOCATION]

        worst <- nPOINTS_SIMPLEX

        ## generate new point for simplex

        ## first extrapolate by a factor -1 through the face of the simplex
        ## across from the high point,i.e.,reflect the simplex from the high point
        parRef <- simplexStep(SIMPLEX, FAC = -1)
        fitRef <- costFunction(FUN, parRef, ...)

        ## check the result
        if (fitRef <= SIMPLEX_FITNESS[1]) {
          ## gives a result better than the best point,so try an additional
          ## extrapolation by a factor 2
          parRefEx <- simplexStep(SIMPLEX, FAC = -2)
          fitRefEx <- costFunction(FUN, parRefEx, ...)
          if (fitRefEx < fitRef) {
            SIMPLEX[worst, ] <- parRefEx
            SIMPLEX_FITNESS[worst] <- fitRefEx
            ALGOSTEP <- "reflection and expansion"
          } else {
            SIMPLEX[worst, ] <- parRef
            SIMPLEX_FITNESS[worst] <- fitRef
            ALGOSTEP <- "reflection"
          }
        } else if (fitRef >= SIMPLEX_FITNESS[worst - 1]) {
          ## the reflected point is worse than the second-highest, so look
          ## for an intermediate lower point, i.e., do a one-dimensional
          ## contraction
          parCon <- simplexStep(SIMPLEX, FAC = -0.5)
          fitCon <- costFunction(FUN, parCon, ...)
          if (fitCon < SIMPLEX_FITNESS[worst]) {
            SIMPLEX[worst, ] <- parCon
            SIMPLEX_FITNESS[worst] <- fitCon
            ALGOSTEP <- "one dimensional contraction"
          } else {
            ## can't seem to get rid of that high point, so better contract
            ## around the lowest (best) point
            SIMPLEX <- (SIMPLEX + rep(SIMPLEX[1, ], each = nPOINTS_SIMPLEX)) / 2
            for (k in 2:NDIM) {
              SIMPLEX_FITNESS[k] <- costFunction(FUN, SIMPLEX[k, ], ...)
            }
            ALGOSTEP <- "multiple contraction"
          }
        } else {
          ## if better than second-highest point, use this point
          SIMPLEX[worst, ] <- parRef
          SIMPLEX_FITNESS[worst] <- fitRef
          ALGOSTEP <- "reflection"
        }

        if (trace >= 3) {
          message(ALGOSTEP)
        }

        ## replace the simplex into the complex
        COMPLEX[LOCATION, ] <- SIMPLEX
        COMPLEX_FITNESS[LOCATION] <- SIMPLEX_FITNESS

        ## sort the complex
        idx <- order(COMPLEX_FITNESS)
        COMPLEX_FITNESS <- COMPLEX_FITNESS[idx]
        COMPLEX <- COMPLEX[idx, , drop = FALSE]
      }

      ## replace the complex back into the population
      POPULATION[k2, ] <- COMPLEX
      POP.FITNESS[k2] <- COMPLEX_FITNESS
    }

    ## At this point, the population was divided in several complexes, each of which
    ## underwent a number of iteration of the simplex (Metropolis) algorithm. Now,
    ## the points in the population are sorted, the termination criteria are checked
    ## and output is given on the screen if requested.

    ## sort the population
    idx <- order(POP.FITNESS)
    POP.FITNESS <- POP.FITNESS[idx]
    POPULATION <- POPULATION[idx, , drop = FALSE]
    if (returnpop) {
      POP.ALL[, , i] <- POPULATION
    }
    POP.FIT.ALL[i, ] <- POP.FITNESS
    BESTMEM.ALL[i, ] <- POPULATION[1, ]

    curBest <- POP.FITNESS[1]

    ## end the optimization if one of the stopping criteria is met

    prevBestVals <- c(curBest, head(prevBestVals, -1))
    reltol <- control$reltol
    if (all(abs(diff(prevBestVals)) <= reltol * (abs(curBest) + reltol))) {
      EXITMSG <- "Change in solution over [tolsteps] less than specified tolerance (reltol)."
      EXITFLAG <- 0
    }

    ## give user feedback on screen if requested
    if (trace >= 1) {
      if (i == 1) {
        message(" Nr Iter  Nr Fun Eval    Current best function    Current worst function")
      }
      if ((i %% control$REPORT == 0) || (!is.na(EXITFLAG))) {
        message(sprintf(
          " %5.0f     %5.0f             %12.6g              %12.6g",
          i, funevals, min(POP.FITNESS), max(POP.FITNESS)
        ))
        if (trace >= 2) {
          message("parameters: ", toString(signif(POPULATION[1, ], 3)))
        }
      }
    }

    if (!is.na(EXITFLAG)) {
      break
    }

    if ((i >= control$maxit) || (funevals >= control$maxeval)) {
      EXITMSG <- "Maximum number of function evaluations or iterations reached."
      EXITFLAG <- 1
      break
    }

    toc <- as.numeric(Sys.time()) - tic
    if (toc > control$maxtime) {
      EXITMSG <- "Exceeded maximum time."
      EXITFLAG <- 2
      break
    }

    ## go to next iteration
    POP.PREV <- POPULATION
    POP.FIT.PREV <- POP.FITNESS
  }
  if (trace >= 1) {
    message(EXITMSG)
  }

  ## return solution
  obj$par <- POPULATION[1, ]
  obj$value <- POP.FITNESS[1]
  obj$convergence <- EXITFLAG
  obj$message <- EXITMSG

  ## store number of function evaluations
  obj$counts <- funevals
  ## store number of iterations
  obj$iterations <- i
  ## store the amount of time taken
  obj$time <- toc

  if (returnpop) {
    ## store information on the population at each iteration
    obj$POP.ALL <- POP.ALL[, , 1:i]
    dimnames(obj$POP.ALL)[[3]] <- paste("iteration", 1:i)
  }
  obj$POP.FIT.ALL <- POP.FIT.ALL[1:i, ]
  obj$BESTMEM.ALL <- BESTMEM.ALL[1:i, ]

  obj
}
