#' Standard methods for Hydromad model objects
#' 
#' A \code{hydromad} object represents a model, which may be fully specified
#' (calibrated) or be defined only by parameter ranges.  The model
#' specification and parameter values are stored along with the observed input
#' and output time-series data.
#' 
#' Several standard methods are available for \code{hydromad} objects:
#' 
#' (note: these are links to the generic functions only)
#' 
#' \code{\link{update}}, \code{\link{predict}}, \code{\link{fitted}},
#' \code{\link{observed}}, \code{\link{residuals}}, \code{\link{coef}},
#' \code{\link{vcov}}, etc.
#' 
#' The \code{\link[=summary.hydromad]{summary}} and
#' \code{\link[=predict.hydromad]{predict}} methods are documented on different
#' pages.
#' 
#' The main plot methods are \code{\link{xyplot.hydromad}} and
#' \code{\link{qqmath.hydromad}}.
#' 
#' \code{isValidModel()} returns \code{TRUE} only if the supplied
#' \code{hydromad} object is fully specified and has a calculated output
#' series.
#' 
#' To help sample parameters, \code{getFreeParsRanges} returns the list of
#' ranges of parameters that do not have fixed parameter values. In particular,
#' it is used in conjunction with \code{\link{evalPars}} to perform sensitivity
#' analyses.
#' 
#' 
#' @aliases hydromad.object update.hydromad fitted.hydromad observed.hydromad
#' residuals.hydromad coef.hydromad getFreeParsRanges vcov.hydromad
#' isValidModel print.hydromad
#' @param object an object of class \code{hydromad}.
#' @param \dots In the \code{update} method, parameter values or ranges for the
#' SMA and/or routing simulation functions can be given, as with the
#' \code{hydromad()} function.
#' @param newdata a \code{\link{ts}}-like object containing a new time series
#' dataset (replacing the original \code{DATA} argument given to the
#' \code{hydromad} function).
#' @param newpars a named list or vector of parameter values; this is
#' equivalent to specifying the same values as named arguments (as in
#' \dQuote{\dots{}}).
#' @param sma,routing,rfit,warmup same arguments as for the
#' \code{\link{hydromad}} function. The \code{update} method allows these to be
#' changed on an existing model object.
#' @param feasible.set,feasible.scores,glue.quantiles the \emph{feasible set}
#' of parameter sets can be specified as a matrix, where parameter values are
#' given in named columns. The corresponding objective function values for each
#' row can be given as \code{feasible.scores}. If \code{glue.quantiles} is
#' omitted or NULL, then overall bounds of the ensemble simulation will be
#' calculated. Otherwise GLUE-like quantiles can be given as
#' \code{glue.quantiles}. See \code{\link{defineFeasibleSet}}.
#' @param and.rescale set to \code{FALSE} to suppress any automatic adjustment
#' of parameters for mass balance.
#' @param which selects either the SMA or routing model, or both models (the
#' default).
#' @param all if \code{TRUE}, return the entire time series for which data
#' exists. Otherwise, the warmup period (specified as an argument to
#' \code{\link{hydromad}} or \code{update}) is stripped off.
#' @param incl.other.vars if \code{TRUE} and model returns a multivariate
#' object (e.g. because of \code{return_components} or \code{return_state}),
#' then return time series for all variables. Otherwise, return only the column
#' named \code{X} or \code{U} (if \code{U=TRUE}). Raises an error if the column
#' is missing.
#' @param feasible.bounds if \code{TRUE}, then ensemble simulation bounds are
#' extracted and returned. This only works if a \emph{feasible set} has been
#' specified using \code{\link{defineFeasibleSet}} or the \code{update} method.
#' Note that the meaning depends on what value of \code{glue.quantiles} was
#' specified to those methods: it might be the overall simulation bounds, or
#' some GLUE-like quantile values. This will be indicated by the returned
#' column names.
#' @param U to return modelled effective rainfall (the output from SMA) rather
#' than streamflow.
#' @param select data series to extract (from the original \code{DATA}
#' argument). Use \code{TRUE} to extract all columns.
#' @param warn by default, \code{coef} gives a warning if the model parameters
#' are not fully specifed (i.e. some parameters have ranges rather than
#' specific values), because it returns a \code{list} rather than a
#' \code{vector} in this case. Setting \code{warn = FALSE} skips the warning.
#' @param etc by default, \code{coef} returns only the model \emph{parameters},
#' which are defined as being numeric and not wrapped in \code{\link{I}()}. If
#' \code{etc = TRUE} is given, then all arguments for the simulation
#' function(s) will be returned, which may include other data types like
#' logicals or time series. In this case the return value is always a
#' \code{list}.
# @param boxcox Placeholder
# @param start Placeholder
#' @return
#' 
#' \code{update} returns a new \code{hydromad} object.
#' 
#' \code{fitted}, \code{observed} and \code{residuals} returns time series.
#' 
#' \code{coef} returns a named numeric vector, or a named \code{list} if one or
#' more parameters are not fully specified.
#' 
#' \code{getFreeParsRanges} returns a named list of parameter ranges, for each
#' parameter that has a range defined. Note that this excludes fixed parameters
#' and parameters defines as a set of discrete values, for which
#' \code{coef(object,warn=FALSE)} should be used instead.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}}, \code{\link{summary.hydromad}},
#' \code{\link{predict.hydromad}}, \code{\link{xyplot.hydromad}},
#' \code{\link{runlist}}, \code{\link{evalPars}}
#' @keywords methods
#' 
#' 
#' 