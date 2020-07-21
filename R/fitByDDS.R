#' Fit a hydromad model using DDS (Dynamically Dimensioned Search) algorithm.
#'
#' Fit a hydromad model using DDS (Dynamically Dimensioned Search) algorithm.
#'
#' This function depends on the \code{ppso} package, available from
#' \href{http://www.rforge.net/ppso/}{http://www.rforge.net/ppso/}. For alternative optimisation algorithms,
#' consider \code{\link{fitBySCE}}.
#'
#' @param MODEL a model specification created by \code{\link{hydromad}}. It
#' should not be fully specified, i.e one or more parameters should be defined
#' by \emph{ranges} of values rather than exact values.
#' @param objective objective function to maximise, given as a
#' \code{function(Q, X, ...)}.  See \code{\link{objFunVal}}.
#' @param control settings for the DDS algorithm. These are the arguments to
#' \code{optim_dds}.
#' @param save Optional \code{function(pars,objective value,model)} that will
#' be called for every model evaluation, for example to save every model run.
#' @return the best model from those sampled, according to the given
#' \code{objective} function. Also, these extra elements are inserted:
#' \item{fit.result}{ the result from \code{\link{SCEoptim}}.  }
#' \item{objective}{ the \code{objective} function used.  } \item{funevals}{
#' total number of evaluations of the model simulation function.  }
#' \item{timing}{ timing vector as returned by \code{system.time}.  }
#' @author Joseph Guillaume \email{josephguillaume@@gmail.com}
#' @seealso \code{optim_dds},\code{\link{fitBySCE}}
#' @references Tolson, B. A., and C. A. Shoemaker (2007) Dynamically
#' dimensioned search algorithm for computationally efficient watershed model
#' calibration, Water Resour. Res., 43, W01413, doi:10.1029/2005WR004723.
#' http://www.agu.org/journals/wr/wr0701/2005WR004723/
#' @keywords optimization
#' @examples
#'
#' data(Cotter)
#' x <- Cotter[1:1000]
#'
#' ## IHACRES CWI model with power law unit hydrograph
#' modx <- hydromad(x, sma = "cwi", routing = "powuh")
#' modx
#'
#' ## run with cut-down settings (for a speedy example only!)
#' foo <- fitByDDS(modx, control = list(
#'   max_number_function_calls = 100,
#'   logfile = NULL,
#'   projectfile = NULL,
#'   load_projectfile = "no"
#' ))
#'
#' summary(foo)
#'
#' ## return value from DDS:
#' str(foo$fit.result)
#' @export
fitByDDS <- function(MODEL, objective = hydromad.getOption("objective"),
                     control = hydromad.getOption("dds.control"), save = NULL) {
  if (!requireNamespace("ppso")) stop("package ppso is required for fitByDDS")
  start_time <- proc.time()
  objective <- buildCachedObjectiveFun(objective, MODEL)
  parlist <- as.list(coef(MODEL, warn = FALSE))
  isok <- sapply(parlist, function(x) !any(is.na(x)))
  parlist <- parlist[isok]
  isfixed <- (sapply(parlist, length) == 1)
  if (all(isfixed)) {
    warning("all parameters are fixed, so can not fit")
    return(MODEL)
  }
  parlist <- parlist[!isfixed]
  ## TODO: if (isTRUE(hydromad.getOption("trace")))
  lower <- sapply(parlist, min)
  upper <- sapply(parlist, max)
  if (is.null(control$initial_estimates)) {
    control$initial_estimates <- as.matrix(sapply(parlist, mean))
  } else {
    if (is.list(control$initial_estimates) | !is.matrix(control$initial_estimates)) stop("control$initial_estimates should be a named column matrix of mode numeric")
    if (!all(names(parlist) %in% rownames(control$initial_estimates))) {
      stop(sprintf(
        "Not all parameters to be optimised have
initial_estimates: %s",
        paste(setdiff(names(parlist), rownames(control$initial_estimates)), collapse = ",")
      ))
    }
    control$initial_estimates <- control$initial_estimates[names(parlist), , drop = F]
  }


  bestModel <- MODEL
  bestFunVal <- -Inf
  do_dds <- function(pars) {
    thisMod <- update(MODEL, newpars = pars)
    if (!isValidModel(thisMod)) {
      return(NA)
    }
    thisVal <- objFunVal(thisMod, objective = objective)
    if (isTRUE(thisVal > bestFunVal)) {
      bestModel <<- thisMod
      bestFunVal <<- thisVal
    }
    if (!is.null(save)) save(pars, thisVal, thisMod)
    return(-thisVal)
  }
  ans <- do.call(
    ppso::optim_dds,
    modifyList(
      control,
      list(
        objective_function = do_dds,
        number_of_parameters = length(control$initial_estimates),
        parameter_bounds = cbind(lower, upper)
      )
    )
  )
  bestModel$msg <- ans$break_flag
  bestModel$funevals <- ans$function_calls
  bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
  bestModel$objective <- objective
  bestModel$fit.call <- match.call()
  bestModel$fit.result <- ans
  return(bestModel)
}
