## Run model and return a vector (usually a timeseries)
## Scalable in length of ts and number of parameters
## Parallelised with results stored on disk
## tempfile is automatically deleted when exiting from R


#' Calculate an objective function on a rolling time series for a matrix of
#' parameters
#'
#' For each row of a named matrix of parameters, run a model and return a
#' vector. By default, \code{evalParsTS} returns the predicted time series.
#' \code{evalParsRollapply} instead calculates an objective function on rolling
#' windows of the resulting time series, i.e. see how the objective function
#' changes over time. Special facilities are provided to evaluate large sample
#' sizes, including with parallelisation.
#'
#' If timeseries are long, then the results matrix will be large
#' (\code{nrow(par.matrix)} x \code{length.out}). By default the results matrix
#' is therefore stored in a \code{ff} file-backed matrix.
#'
#' Individual model evaluations are generally very fast, so parallelisation is
#' only really worthwhile when large numbers of evaluations are needed.
#' Parallelisation method \code{"clusterApply"} uses multiple R sessions that
#' write to a shared \code{ff} object. It only operates on a single
#' multicore machine. Parallelisation method \code{"foreach"} offers a broader
#' range of options but only works if the final results matrix is small enough
#' to fit in memory.
#'
#' @importFrom stats update
#' @importFrom parallel clusterCall
#'
#' @aliases evalParsTS evalParsRollapply
#' @param par.matrix Named matrix or data.frame of parameter values, with each
#' row corresponding to a model realisation to evaluate
#' @param object an object of class \code{hydromad}.
#' @param fun function that takes a hydromad object and returns a vector, by
#' default the fitted timeseries.
#' @param length.out Length of output vector returned by \code{fun}. If
#' missing, \code{fun} will be run on the first parameter set in
#' \code{par.matrix}.
#' @param \dots Additional arguments to \code{fun}
#' @param parallel If \code{"clusterApply"}, evaluate parameters in parallel
#' using a local cluster. The implementation assumes that the \code{ff}
#' file \code{filehash.name} can be written to simultaneously all cluster
#' workers.
#' @param filehash.name Name of \code{ff} file in which to store
#' results, allowing large samples that do not fit in memory. Defaults to
#' \code{tempfile()}, which is automatically deleted when exiting from R. To
#' store results in memory, set \code{filehash.name=NULL}.
#' @param width integer specifying window width, aligned center. Passed to
#' \code{\link[zoo]{rollapply}}
#' @param objective the objective function or expression, which can refer to Q
#' and X.  See \code{\link{objFunVal.hydromad}}
#' @return Either a matrix or \code{ff} file-backed matrix, with each
#' row being a time series of rolling objective functions for the corresponding
#' row of \code{par.matrix}
#' @note When using \code{ff}, performance may be improved by specifying
#' \code{options(ffcaching='mmeachflush')}.
#' @author Joseph Guillaume
#' @seealso \code{\link{evalPars}} to calculate a single objective function for
#' each parameter set
#' @references Herman, J. D., P. M. Reed, and T. Wagener. 2013. "Time-Varying
#' Sensitivity Analysis Clarifies the Effects of Watershed Model Formulation on
#' Model Behavior." Water Resources Research 49 (3): 1400-1414.
#' doi: \href{http://dx.doi.org/10.1002/wrcr.2012410.1002/wrcr.20124}{here}
#' @keywords models
#' @examples
#'
#' data(Cotter)
#' obs <- Cotter[1:1000]
#'
#' ## Define rainfall-runoff model structure
#' object <- hydromad(obs,
#'   sma = "cwi", routing = "expuh",
#'   tau_q = c(0, 2), tau_s = c(2, 100), v_s = c(0, 1)
#' )
#'
#' ## Set the random seed to obtain replicable results
#' set.seed(19)
#'
#' # Draw 10 Latin Hypercube random samples
#' par.matrix <- parameterSets(getFreeParsRanges(object), samples = 10)
#' # Calculate rolling time series of r.squared for each parameter set,
#' #   keeping results in memory
#' runs <- evalParsRollapply(par.matrix, object,
#'   objective = hmadstat("r.squared"), filehash.name = NULL
#' )
#' \dontrun{
#' ## Setup parallelisation on three cores
#' library(parallel)
#' hydromad.options("parallel" = list("evalParsTS" = "clusterApply"))
#' cl <- makeCluster(3)
#' clusterEvalQ(cl, library(hydromad))
#'
#' par.matrix <- parameterSets(getFreeParsRanges(object), samples = 1000)
#' # Calculate rolling time series of r.squared for each parameter set,
#' #  storing result in tempfile()
#' # Takes about 2 minutes
#' runs <- evalParsRollapply(par.matrix, object,
#'   objective = hmadstat("r.squared")
#' )
#'
#' # Excerpt of results
#' runs
#' # Path of backing file - about 7MB
#' filename(runs)
#' # ff object can be used like a regular matrix, e.g. plotting the
#' #  rolling time series of R.squared for the first parameter set
#' plot(runs[1, ])
#'
#' ## Do the same with foreach
#' library(doParallel)
#' registerDoParallel(cl)
#' hydromad.options("parallel" = list("evalParsTS" = "foreach"))
#' runs <- evalParsRollapply(par.matrix, object,
#'   objective = hmadstat("r.squared"), filehash.name = NULL
#' )
#' ## runs is a matrix
#' }
#'
#' @export
## Calculate objective function on a rolling window
evalParsRollapply <- function(par.matrix, object,
                              width = 30,
                              objective = hydromad.getOption("objective"),
                              parallel = hydromad.getOption("parallel")[["evalParsTS"]],
                              filehash.name = tempfile()) {
  fun <- function(thisMod, width, objective) {
    rollapply(cbind(Q = observed(thisMod), X = fitted(thisMod)),
      width = width, by.column = FALSE,
      FUN = objFunVal, objective = objective
    )
  }

  evalParsTS(par.matrix, object, fun,
    length.out = length(observed(object)) - width + 1,
    parallel = parallel,
    filehash.name = filehash.name,
    width = width, objective = objective
  )
}



#' @rdname evalParsRollapply
#' @export
evalParsTS <- function(par.matrix, object,
                       fun = function(thisMod) fitted(thisMod),
                       length.out = NULL,
                       ...,
                       parallel = hydromad.getOption("parallel")[["evalParsTS"]],
                       filehash.name = tempfile()) {
  .args <- list(...)

  if (is.null(length.out)) {
    warning("length.out missing, running model to evaluate it")
    res <- fun(update(object, newpars = par.matrix[1, ]), ...)
    length.out <- length(res)
    cat("fun returns vector with length.out=", length.out, "\n")
  }

  # Sets default settings for parallelisation if missing
  parallel <- hydromad.parallel(parallel)

  if (parallel$method == "foreach" && !is.null(filehash.name)) {
    filehash.name <- NULL
    warning("ignoring filehash.name, 'foreach' parallelisation does not support writing directly to disk")
  }

  if (is.null(filehash.name)) {
    ## Don't use disk (ff)
    results <- matrix(NA, nrow = nrow(par.matrix), ncol = length.out)
    if (parallel$method == "clusterApply") {
      warning("setting parallel$method='none', 'clusterApply' parallelisation requires filehash.name to be non-null")
      parallel$method <- "none"
    }
  } else {
    ## Use disk
    if (!requireNamespace("parallel")) stop("package parallel is required for evalParsTS if filehash.name is not NULL and parallel$method is not 'foreach'")
    if (!requireNamespace("ff")) stop("package ff is required for evalParsTS if filehash.name is not NULL and parallel$method is not 'foreach'")
    results <- ff::ff(vmode = "double", dim = c(nrow(par.matrix), length.out), filename = filehash.name)
  }

  ## Do runs, storing all ts
  cat(sprintf("Running %d model evaluations with parallelisation='%s'\n", nrow(par.matrix), parallel$method))
  switch(parallel$method,
    "foreach" = {
      opts <- hydromad.options()
      export <- parallel$export
      results <- foreach::foreach(
        p = iterators::iter(par.matrix, by = "row"),
        .packages = parallel$packages,
        .inorder = TRUE,
        .export = export,
        .final = function(x) do.call(rbind, lapply(x, coredata)),
        .options.redis = list(async = parallel$async)
      ) %dopar% {
        # Work-around for hydromad functions to have access to .export
        for (e in export) assign(e, get(e), envir = .GlobalEnv)
        # Work-around to use same opts as in user's environment
        hydromad.options(opts)
        thisMod <- update(object, newpars = p)
        if (!isValidModel(thisMod)) {
          return(NULL)
        }
        do.call(fun, modifyList(.args, list(thisMod = thisMod)))
      }
    },
    "clusterApply" = {
      lapply(c("ff", parallel$packages), function(pkg) parallel::clusterCall(cl, library, pkg, character.only = TRUE))
      if (length(parallel$export) > 0) parallel::clusterExport(cl, parallel$export)
      parallel::clusterExport(cl, c("object"), envir = environment())
      parallel::parLapply(
        cl = cl, as.list(1:nrow(par.matrix)),
        function(ip) {
          thisMod <- update(object, newpars = par.matrix[ip, ])
          if (!isValidModel(thisMod)) {
            return(NULL)
          }
          results[ip, ] <- do.call(fun, modifyList(.args, list(thisMod = thisMod)))
          ip
        }
      )
    },
    {
      lapply(
        as.list(1:nrow(par.matrix)),
        function(ip) {
          thisMod <- update(object, newpars = par.matrix[ip, ])
          if (!isValidModel(thisMod)) {
            return(NULL)
          }
          results[ip, ] <<- do.call(fun, modifyList(.args, list(thisMod = thisMod)))
          ip
        }
      )
    }
  ) ## switch parallel
  ## ff matrix of ts
  results
}

#' @import utils
utils::globalVariables(c("cl", "p"))
