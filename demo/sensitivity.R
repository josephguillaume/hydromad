#' Support for sensitivity analysis in hydromad
#'
#' Sensitivity analysis can be performed on hydromad objects using the
#' \code{\link[sensitivity]{sensitivity}} package. The sensitivity analysis
#' functions (e.g. \code{morris}, \code{sobol2007}) are called directly,
#' specifying the arguments \code{model}, \code{object} and \code{objective} as
#' described below.
#'
#'
#' @aliases hydromad_sensitivity sensitivity
#' @param model Must be \code{model=\link{evalPars}}, which evaluates the
#' \code{objective} using the model specified by \code{object}
#' @param object an object of class \code{hydromad}.
#' @param objective the objective function or expression, which can refer to Q
#' and X.  See \code{\link{objFunVal.hydromad}}
#' @param parallel (optional) parallelisation of model runs, see
#' \code{\link{evalPars}} and \code{\link{hydromad_parallelisation}}
#' @author Joseph Guillaume
#' @seealso Utilities \code{\link{evalPars}} and
#' \code{\link{getFreeParsRanges}}. \code{\link{tellTS}} for calculating
#' sensitivity indices on a time series.
#' @references Shin, Mun-Ju, Joseph H.A. Guillaume, Barry F.W. Croke, and
#' Anthony J. Jakeman. 2013. "Addressing Ten Questions about Conceptual
#' Rainfall-runoff Models with Global Sensitivity Analyses in R." Journal of
#' Hydrology 503 (October): 135-52.
#' doi:\href{http://dx.doi.org/10.1016/j.jhydrol.2013.08.04710.1016/j.jhydrol.2013.08.047}
#' @keywords models
#' @examples
#'
#' library(sensitivity)
#'
#' ## Load data
#' data(Cotter)
#' obs <- Cotter[1:1000]
#'
#' ## Define rainfall-runoff model structure
#' model.str <- hydromad(obs, sma = "cwi", routing = "expuh",
#'    tau_s = c(2,100), v_s = c(0,1))
#'
#'
#' ## Set the random seed to obtain replicable results
#' set.seed(19)
#'
#'
#' ################################################################################
#' ## Sensitivity using Morris method of NSE* objective function to
#' ##  IHACRES-CWI model parameters using a subset of data from Cotter catchment
#'
#' ## Run Morris Sensitivity analysis
#' mm <- morris(
#'              ## Utility function to evaluate parameters on a model
#'              model=evalPars,
#'              ## Names of factors/parameters
#'              factors=names(getFreeParsRanges(model.str)),
#'              ## Number of repetitions
#'              r=4,
#'              ## List specifying design type and its parameters
#'              design=list(type="oat",levels=10,grid.jump=2),
#'              ## Minimum value of each non-fixed parameter
#'              binf=sapply(getFreeParsRanges(model.str),min),
#'              ## Maximum value of each non-fixed parameter
#'              bsup=sapply(getFreeParsRanges(model.str),max),
#'              ## Hydromad model object to use to evaluate parameters
#'              object=model.str,
#'              ## NSE* objective function
#'              objective=~ hmadstat("r.squared")(Q, X) /
#'                 (2-hmadstat("r.squared")(Q, X))
#'              )
#'
#' print(mm)
#'
#' ## Default plot of mu.star (mean absolute sensitivity) and sigma (sd of sensitivity)
#' plot(mm,main = "Sensitivity NSE*~IHACRES-CWI parameters with Cotter data")
#'
#' ## For custom plots, mu.star and sigma can be explicitly calculated
#' mu.star <- apply(mm$ee, 2, function(x) mean(abs(x)))
#' sigma <- apply(mm$ee, 2, sd)
#' plot(mu.star, sigma, pch = 20, xlab = expression(mu^"*"),
#'      ylab = expression(sigma))
#' text(mu.star, sigma, labels = colnames(mm$ee), pos = 4,offset=0.4)
#'
#'
#' #' ################################################################################
#' ## We can perform sensitivity analysis on the prediction function, to identify
#' ## whether parameters that are insensitive with the objective function can be fixed:
#' ##  Sensitivity using Morris method of F20 prediction function
#' ##  (frequency/ratio of days of flow below 20%ile)
#' ##  to IHACRES-CWI model parameters using a subset of data from Cotter catchment
#'
#' ## Calculate observed 20%ile flow
#' thres.Q20 <- as.numeric(quantile(obs$Q, probs = c(0.2), na.rm = TRUE))
#' ## Define function to calculate number of days below thres.Q20
#' F20 <- function(X) length(which(X < thres.Q20)) / length(X)
#'
#' mm <- morris(
#'   model = evalPars,
#'   factors = names(getFreeParsRanges(modx)),
#'   r = 4,
#'   design = list(
#'     type = "oat",
#'     levels = 10,
#'     grid.jump = 2
#'   ),
#'   binf = sapply(getFreeParsRanges(modx), min),
#'   bsup = sapply(getFreeParsRanges(modx), max),
#'   object = modx,
#'  ## Change objective to use the prediction function defined above
#'  ## See ?objFunVal
#'  ## The objective function can refer to Q and X,
#'  ##  representing observed and modelled flow, respectively.
#'  ##  It should return a single numeric value.
#'   objective = ~ F20(X)
#' )
#'
#' print(mm)
#' plot(mm, main = "Sensitivity F20~IHACRES-CWI parameters with Cotter data")
#'
#' ## See the Morris (1991) journal paper(?morris) for a more comprehensive interpretation.
#' ## These results are only an example, but if sampling were sufficient, the results would
#' ##  suggest that NSE* is insensitive to tau_s (mu.star=0.065), but F20 is sensitive to
#' ##  tau_s (mu.star=0.34)
#' ## The parameter value selected therefore affects the prediction, so tau_s cannot be
#' ##  fixed to improve identifiability unless its value is sufficiently certain.
#' ## Using a different objective function, e.g. NSElog* might improve identifiability.
#' ## The settings for Morris used may also influence results. Try changing: the parameter ranges,
#' ##  number of replicates, number of levels and length of data series.
#' ## A SOBOL sensitivity analysis would also provide confidence intervals around the TSI.

#'
#' ################################################################################
#'
#' ## Sensitivity using SOBOL2002 method of NSElog* prediction function to
#' ##  IHACRES-CWI model parameters using a subset of data from Cotter catchment
#' ## This might take a while, potentially ~15min
#'
#' \dontrun{
#' n <- 1000 ## Number of random samples for parameters
#' ss <- sobol2002(model = evalPars,
#'                 ## Draw two random samples of parameters
#'                 X1 = parameterSets(getFreeParsRanges(model.str),n),
#'                 X2 = parameterSets(getFreeParsRanges(model.str),n),
#'                 ## Number of bootstrap replicates
#'                 nboot = 100,
#'                 ## Hydromad model object to use to evaluate parameters
#'                 object=model.str,
#'                 ## NSElog* objective function (using the logarithm of Q and X)
#'                 objective=~ hmadstat("r.sq.log")(Q, X) /
#'                   (2-hmadstat("r.sq.log")(Q, X))
#'                 )
#'
#' ## Show results
#' print(ss)
#' plot(ss)
#' }
#'
#'
#' ## tau_s appears to still be insensitive with NSElog* (TSI=0.127), suggesting
#' ## the information in this dataset may not be sufficient to identify this parameter.
#'
