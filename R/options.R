## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

.defaultHydromadOptions <- function() {
  list(
    sma = NULL,
    routing = NULL,


#' Simple constant runoff proportion
#' 
#' Simple constant runoff proportion: a constant fraction of rainfall reaches
#' the stream.
#' 
#' 
#' @aliases scalar scalar.sim absorbScale.hydromad.scalar
#' @param DATA time-series-like object with a column P (precipitation).
#' @param scale fraction of rainfall that becomes effective.  If this parameter
#' is set to \code{NA} (as it is by default) in \code{\link{hydromad}} it will
#' be set by mass balance calculation.
#' @param return_state ignored.
#' @return the simulated effective rainfall, a time series of the same length
#' as the input series.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "scalar")} to work with models as
#' objects (recommended).
#' @keywords models
#' @examples
#' 
#' ## view default parameter ranges:
#' str(hydromad.options("scalar"))
#' 
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "scalar", routing = "expuh")
#' mod0
#' 
#' ## simulate with some arbitrary parameter values
#' testQ <- predict(update(mod0, scale = 0.5, tau_s = 10))
#' xyplot(cbind(HydroTestData[,1:2], scalar.Q = testQ))
#' 
    scalar = scalar.ranges(),
    cwi = cwi.ranges(),
    cmd = cmd.ranges(),


#' Sacramento Soil Moisture Accounting model
#' 
#' Sacramento Soil Moisture Accounting model.  Developed by the US National
#' Weather Service.
#' 
#' This description of the model is given by Burnash (1995):
#' 
#' \dQuote{The moisture accounting system utilized in the Sacramento Catchment
#' Model is a carefully structured representation of the catchment's soil
#' moisture storage system. It is based on using simple approximations of many
#' of those soil moisture processes which have been reported in the hydrologic
#' literature. The authors have organised these approximations in a manner
#' which would allow the determination of many catchment characteristics from
#' carefully selected portions of the catchment's hydrologic record. Inasmuch
#' as many of the catchment characteristics are related to the soil moisture
#' capabilities of the catchment, an intelligent application of the model start
#' with a good understanding of the three basic types of soil moisture which
#' can potentially influence catchment runoff conditions. These soil moisture
#' types are: (1) Hygroscopic Water, (2) Tension Water and (3) Free Water. }
#' 
#' [...]
#' 
#' \dQuote{Streamflow as computed by the Sacramento Catchment Model is the
#' result of processing precipiatation through an algorithm representing the
#' uppermost soil mantle identified as the upper zone and a deeper portion of
#' the soil mantle or lower zone. The algorithm computes runoff in five basic
#' forms. These are (1) direct runoff from permanant and temporary impervious
#' areas, (2) surface runoff due to precipitation occurring at a rate faster
#' than percolation and interflow can take place when both upper zone storages
#' are full, (3) interflow resulting from the lateral drainage of a temporary
#' free water storage, (4) supplemental base flow, and (5) primary base flow.}
#' (Burnash, 1995)
#' 
#' The default parameter ranges were taken from Blasone et. al. (2008).
#' 
#' Note that the Sacramento model potentially suffers from numerical
#' instabilities, which can be seen for example as discontinuities in output
#' and derivatives of outputs (see Hendrickson et al. 1988). Ideally, the
#' underlying differential equations of the model would be solved using a
#' numerically robust timestepping scheme (see Clark & Kavetski 2010). The
#' hydromad package makes use of an existing implementation. To help remedy the
#' numerical instability, the argument \code{min_ninc} has been added, which
#' defines the minimum number of inner loops used within each timestep. The
#' user is encouraged to test the effect of increasing \code{min_ninc} on their
#' dataset.
#' 
#' @aliases sacramento sacramento.sim
#' @param DATA time-series-like object with columns \code{P} (precipitation,
#' mm) and \code{E} (potential evapo-transpiration, mm, scaled by
#' \code{etmult}).
#' @param uztwm Upper zone tension water maximum capacity (mm).
#' @param uzfwm Upper zone free water maximum capacity (mm).
#' @param uzk Lateral drainage rate of upper zone free water expressed as a
#' fraction of contents per day.
#' @param pctim The fraction of the catchment which produces impervious runoff
#' during low flow conditions.
#' @param adimp The additional fraction of the catchment which exhibits
#' impervious characteristics when the catchment's tension water requirements
#' are met.
#' @param zperc Maximum percolation (from upper zone free water into the lower
#' zone) rate coefficient.
#' @param rexp An exponent determining the rate of change of the percolation
#' rate with changing lower zone water contents.
#' @param lztwm Lower zone tension water maximum capacity (mm).
#' @param lzfsm Lower zone supplemental free water maximum capacity (mm).
#' @param lzfpm Lower zone primary free water maximum capacity (mm).
#' @param lzsk Lateral drainage rate of lower zone supplemental free water
#' expressed as a fraction of contents per day.
#' @param lzpk Lateral drainage rate of lower zone primary free water expressed
#' as a fraction of contents per day.
#' @param pfree Direct percolation fraction from upper to lower zone free water
#' (the percentage of percolated water which is available to the lower zone
#' free water aquifers before all lower zone tension water deficiencies are
#' satisfied).
#' @param etmult Multiplier applied to \code{DATA$E} to estimate potential
#' evapotranspiration.
#' @param dt Length of each time step in days.
#' @param uztwc_0 Initial upper zone tension water contents as proportion of
#' \code{uztwm}
#' @param uzfwc_0 Initial upper zone free water content as proportion of
#' \code{uzfwm}
#' @param lztwc_0 Initial lower zone tension water content as proportion of
#' \code{lztwm}
#' @param lzfsc_0 Initial lower zone free water secondary as proportion of
#' \code{lzfsm}
#' @param lzfpc_0 Initial lower zone free water primary as proportion of
#' \code{lzfpm}
#' @param adimc_0 Initial additional impervious flow store, as proportion of
#' \code{uztwm+lztwm}
#' @param min_ninc Minimum number of inner iterations. This is a simple attempt
#' to improve numerical stability. See Details.
#' @param return_state to return time series of each state variable and flow
#' component
#' @return the simulated effective rainfall (\dQuote{total channel inflow}), a
#' time series of the same length as the input series.
#' 
#' if \code{return_state=TRUE}, a list with components: \item{uztwc}{Upper zone
#' tension water content} \item{uzfwc}{Upper zone free water content}
#' \item{lztwc}{Lower zone tension water content} \item{lzfsc}{Lower zone free
#' secondary water content} \item{lzfpc}{Lower zone free primary water content}
#' \item{adimc}{Tension water contents of the additional impervious area}
#' \item{sett}{Cumulative total evapotranspiration} \item{se1}{Cumulative
#' evapotranspiration from upper zone tension water} \item{se3}{Cumulative
#' evapotranspiration from lower zone tension water} \item{se4}{Cumulative
#' evapotranspiration} \item{se5}{Cumulative evapotranspiration from riparian
#' zone} \item{roimp}{Runoff from impervious area} \item{sdro}{Six hour sum of
#' runoff (?)} \item{ssur}{Surface runoff} \item{sif}{Interflow}
#' \item{bfp}{Primary baseflow} \item{bfs}{Secondary baseflow}
#' \item{bfcc}{Channel baseflow (bfp+bfs)}
#' @author Felix Andrews \email{felix@@nfrac.org} and Joseph Guillaume, based
#' on code from the University of Arizona MOSCEM project
#' @seealso \code{\link{hydromad}(sma = "sacramento")} to work with models as
#' objects (recommended).
#' @references Burnash, R.J.C (1995). The NWS River Forecast System --
#' Catchment Modeling.  In: Vijay P. Singh (ed.), \emph{Computer models of
#' watershed hydrology.} Revised edition, Highlands Ranch, Colo. : Water
#' Resources Publications, c1995.  \url{http://www.wrpllc.com/books/cmwh.html}.
#' 
#' Blasone, R., J.A. Vrugt, H. Madsen, D. Rosbjerg, B.A. Robinson, G.A.
#' Zyvoloski (2008). Generalized likelihood uncertainty estimation (GLUE) using
#' adaptive Markov Chain Monte Carlo sampling. \emph{Advances in Water
#' Resources} 31, pp. 630-648.
#' 
#' Hendrickson, Jene' D., Soroosh Sorooshian, and Larry E. Brazil (1988)
#' Comparison of Newton-Type and Direct Search Algorithms for Calibration of
#' Conceptual Rainfall-Runoff Models. \emph{Water Resources Research} 24 (5):
#' 691-700.  \url{http://dx.doi.org/10.1029/WR024i005p00691}
#' 
#' Clark, Martyn P., and Dmitri Kavetski (2010) Ancient Numerical Daemons of
#' Conceptual Hydrological Modeling: 1. Fidelity and Efficiency of Time
#' Stepping Schemes.” Water Resources Research 46 (10).
#' \url{http://dx.doi.org/10.1029/2009WR008894}
#' @keywords models
#' @examples
#' 
#' ## view default parameter ranges:
#' str(hydromad.options("sacramento"))
#' 
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "sacramento")
#' mod0
#' 
#' ## simulate with some arbitrary parameter values
#' set.seed(2)
#' mod1 <- simulate(update(mod0, etmult = 0.01), 1, sampletype =
#' "random")[[1]]
#' 
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(window(cbind(HydroTestData[,1:2], sacramento = testQ), start = 100))
#' mod1
#' 
#' ## show effect of increase/decrease in each parameter
#' parRanges <- hydromad.getOption("sacramento")
#' parsims <- mapply(val = parRanges, nm = names(parRanges),
#'   FUN = function(val, nm) {
#'     lopar <- min(val)
#'     hipar <- max(val)
#'     names(lopar) <- names(hipar) <- nm
#'     fitted(runlist(decrease = update(mod1, newpars = lopar),
#'                    increase = update(mod1, newpars = hipar)))
#'   }, SIMPLIFY = FALSE)
#' 
#' xyplot.list(parsims, superpose = TRUE, layout = c(1,NA),
#'             strip = FALSE, strip.left = TRUE,
#'             main = "Simple parameter perturbation example") +
#'   latticeExtra::layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
#' 
    sacramento = sacramento.ranges(),


#' Single-bucket Soil Moisture Accounting models
#' 
#' Single-bucket Soil Moisture Accounting models with saturated/unsaturated
#' zones and interception.
#' 
#' From formulations given in Bai et. al. (2009), which were based on Farmer
#' et. al. (2003).
#' 
#' The general mass balance structure is: \deqn{dS/dt = p - q(S) - e(S, Ep)}
#' 
#' The default parameter ranges were also taken from Bai et. al. (2009).
#' 
#' @aliases bucket bucket.sim
#' @param DATA time-series-like object with columns P (precipitation, mm) and E
#' (potential evapo-transpiration, mm).
#' @param Sb Maximum soil water storage (mm).
#' @param fc Field capacity (0 - 1).
#' @param a.ei Interception coefficient (\eqn{\alpha_{ei}}).
#' @param M Fraction of catchment area covered by deep rooted vegetation.
#' @param a.ss Recession coefficients for subsurface flow from saturated zone
#' (\eqn{\alpha_{ss}}).
#' @param etmult Multiplier for the \code{E} input data.
#' @param S_0 Initial soil moisture level as a fraction of \code{Sb}.
#' @param return_state to return the series U, S and ET (evapotranspiration).
#' @return the simulated effective rainfall, a time series of the same length
#' as the input series.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "bucket")} to work with models as
#' objects (recommended).
#' @references Farmer, D., M. Sivapalan, Farmer, D. (2003). Climate, soil and
#' vegetation controls upon the variability of water balance in temperate and
#' semiarid landscapes: downward approach to water balance analysis.
#' \emph{Water Resources Research} 39(2), p 1035.
#' 
#' Bai, Y., T. Wagener, P. Reed (2009). A top-down framework for watershed
#' model evaluation and selection under uncertainty. \emph{Environmental
#' Modelling and Software} 24(8), pp. 901-916.
#' @keywords models
#' @examples
#' 
#' ## view default parameter ranges:
#' str(hydromad.options("bucket"))
#' 
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "bucket", routing = "expuh")
#' mod0
#' 
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, Sb = 10, fc = 0.5, M = 0.5, etmult = 0.05,
#'                      a.ei = 0.05, a.ss = 0.01, tau_s = 10)
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[,1:2], bucket = testQ))
#' 
#' ## show effect of increase/decrease in each parameter
#' parRanges <- hydromad.getOption("bucket")
#' parsims <- mapply(val = parRanges, nm = names(parRanges),
#'   FUN = function(val, nm) {
#'     lopar <- min(val)
#'     hipar <- max(val)
#'     names(lopar) <- names(hipar) <- nm
#'     fitted(runlist(decrease = update(mod1, newpars = lopar),
#'                    increase = update(mod1, newpars = hipar)))
#'   }, SIMPLIFY = FALSE)
#' 
#' xyplot.list(parsims, superpose = TRUE, layout = c(1,NA),
#'             strip = FALSE, strip.left = TRUE,
#'             main = "Simple parameter perturbation example") +
#'   layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
#' 
    bucket = bucket.ranges(),


#' GR4J rainfall runoff model
#' 
#' GR4J model (mode`le du Ge´nie Rural a` 4 parame`tres Journalier).
#' 
#' The default parameter ranges were taken from the "80% confidence intervals"
#' given in Perrin et. al. (2003).
#' 
#' @aliases gr4j gr4j.sim gr4jrouting gr4jrouting.sim
#' @param DATA time-series-like object with columns P (precipitation, mm) and E
#' (potential evapo-transpiration, mm).
#' @param U effective rainfall series.
#' @param x1 maximum capacity of the production store (mm).
#' @param x2 groundwater exchange coefficient (mm).
#' @param x3 one day ahead maximum capacity of the routing store (mm).
#' @param x4 time base of unit hydrograph UH1 (time steps).
#' @param etmult Multiplier for the \code{E} input data.
#' @param S_0 Initial soil moisture level as fraction of \code{x1}.
#' @param R_0 Initial groundwater reservoir level as fraction of \code{x3}.
#' @param split Fraction to go into quick flow routing, usually fixed at 0.9.
#' @param return_state to return the series U, S (storage) and ET
#' (evapotranspiration).
#' @param return_components to return the series Xr, Xd and R (reservoir
#' level).
#' @param epsilon values smaller than this in the output will be set to zero.
#' @param transformed transform parameters before use to improve
#' identifiability. They can be untransformed using
#' \code{\link{gr4j.transformpar}}
#' @return the simulated effective rainfall, a time series of the same length
#' as the input series.
#' @author Felix Andrews \email{felix@@nfrac.org} and Joseph Guillaume
#' \email{josephguillaume@@gmail.com}
#' @seealso \code{\link{hydromad}(sma = "gr4j", routing = "gr4jrouting")} to
#' work with models as objects (recommended).
#' @references Perrin, C., C. Michel, et al. (2003). "Improvement of a
#' parsimonious model for streamflow simulation." \emph{Journal of Hydrology}
#' 279(1-4): 275-289
#' 
#' \url{http://www.cemagref.fr/webgr/Modelesgb/gr4j/fonctionnement_gr4jgb.htm}
#' @keywords models
#' @examples
#' 
#' ## view default parameter ranges:
#' str(c(hydromad.getOption("gr4j"),
#'       hydromad.getOption("gr4jrouting")))
#' 
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "gr4j", routing = "gr4jrouting")
#' mod0
#' 
#' ## example from
#' ## http://www.cemagref.fr/webgr/Scilab/CONT_EN/HELP_HYDROGR/c_GR4J.htm
#' dat <-
#'   cbind(P = c(0,0,0,0,0,0.04,0.59,0.03,0.01,0.16,0.37,8.76,2.65,
#'           0.05,0.02,0.02,0.38,0.00,0.02,0.46,4.46,7.71,5.71,0.79,1.33),
#'         E = c(0,0,0,0,0,0.24,0.24,0.24,0.24,0.24,0.25,0.25,0.26,
#'           0.27,0.28,0.32,0.33,0.34,0.35,0.36,0.36,0.37,0.37,0.38,0.38))
#' datz <- zoo(dat, as.Date("2000-01-01") + 1:nrow(dat))
#' modz <- hydromad(datz, sma = "gr4j", routing = "gr4jrouting",
#'     x1 = 665, x2 = 1.18, x3 = 90, x4 = 3.8, S_0 = 0.6, R_0 = 0.7)
#' xyplot(predict(modz, return_state = TRUE, return_components = TRUE),
#'        strip = FALSE, strip.left = TRUE)
#' 
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, x1 = 100, x2 = 20, x3 = 1, x4 = 10)
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[,1:2], gr4j = testQ))
#' 
#' ############################################
#' ## show effect of increase/decrease in each parameter
#' parRanges <- c(hydromad.getOption("gr4j")[1],
#'                hydromad.getOption("gr4jrouting"))
#' parsims <- mapply(val = parRanges, nm = names(parRanges),
#'   FUN = function(val, nm) {
#'     lopar <- min(val)
#'     hipar <- max(val)
#'     names(lopar) <- names(hipar) <- nm
#'     fitted(runlist(decrease = update(mod1, newpars = lopar),
#'                    increase = update(mod1, newpars = hipar)))
#'   }, SIMPLIFY = FALSE)
#' 
#' xyplot.list(parsims, superpose = TRUE, layout = c(1,NA),
#'             strip = FALSE, strip.left = TRUE,
#'             main = "Simple parameter perturbation example") +
#'   layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
#' 
#' 
#' ############################################
#' # Example optimisation, using transformed parameters
#' 
#' data(Cotter)
#' x <- Cotter[1:1000]
#' 
#' #Specify gr4j model
#' mod0 <- hydromad(x, sma = "gr4j", routing = "gr4jrouting",transformed=TRUE)
#' #Use transformed parameter ranges
#' mod0 <- update(mod0,newpars=gr4j.transformpar(c(hydromad.getOption("gr4j"),
#'                                                 hydromad.getOption("gr4jrouting")
#'                                                 )))
#' #Allow etmult to vary, because we're using temperature data instead of PET.
#' mod0<-update(mod0,etmult=c(0.05,1.5))
#' # Broaden a single parameter range, just as an example
#' mod0<-update(mod0,x1=gr4j.transformpar(list(x1=c(100,5000)))[["x1"]])
#'                
#' mod0
#' 
#' ## now try to fit the free parameters
#' set.seed(10)
#' fit1 <- fitByOptim(mod0)
#' 
#' fit1
#' summary(fit1)
#' xyplot(fit1)
#' 
#' #Parameters in original parameter space
#' gr4j.transformpar(coef(fit1),back=T)
#' 
#' 
#' 
    gr4j = gr4j.ranges(),
    gr4jrouting = gr4jrouting.ranges(),


#' Australian Water Balance Model (AWBM)
#' 
#' Australian Water Balance Model (AWBM): simple 3 bucket model.
#' 
#' This is a very simple model: saturation excess from three buckets with
#' different capacities are added together with fractional areas for weights.
#' 
#' This is the soil moisture accounting component; the model described in the
#' literature has a two-store routing component also, with parameters
#' \var{BFI}, \eqn{K_b} and \eqn{K_s}. These correspond directly to the
#' \code{\link{expuh}} routing model parameters \code{v_s}, \code{tau_s} and
#' \code{tau_q}, so the full model can be specified as:
#' 
#' \code{hydromad(..., sma = "awbm", routing = "expuh", rfit = list("sriv",
#' order = c(2, 1)))}
#' 
#' @aliases awbm awbm.sim
#' @param DATA time-series-like object with columns P (precipitation, mm) and E
#' (potential evapo-transpiration, mm).
#' @param cap.ave average soil water storage capacity (mm). This is not used
#' directly in the model, but rather to define default values of the other
#' parameters.
#' @param cap1,cap2,cap3 soil water storage capacities (mm) for the three
#' fractional areas.
#' @param area1,area2,area3 fractional areas with corresponding capacities.
#' These must sum to 1.
#' @param etmult multiplier for the \code{E} input data.
#' @param S1_0,S2_0,S3_0 initial soil moisture levels (mm).
#' @param return_state to return the soil moisture levels \code{S1}, \code{S2}
#' and \code{S3} together with effective rainfall \code{U}.
#' @return the simulated effective rainfall, a time series of the same length
#' as the input series.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "awbm")} to work with models as objects
#' (recommended).
#' @references Boughton, W. (2004). The australian water balance model.
#' Environmental Modelling & Software 19 (10), 943-956.
#' \url{http://dx.doi.org/10.1016/j.envsoft.2003.10.007}
#' @keywords models
#' @examples
#' 
#' ## view default parameter ranges:
#' str(hydromad.options("awbm"))
#' 
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "awbm", routing = "expuh")
#' mod0
#' 
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, cap.ave = 40, etmult = 0.1, tau_s = 10)
#' 
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[,1:2], awbm = testQ))
#' 
#' ## show effect of increase/decrease in each parameter
#' parRanges <- list(cap.ave = c(1, 1000), etmult = c(0.01, 1))
#' parsims <- mapply(val = parRanges, nm = names(parRanges),
#'   FUN = function(val, nm) {
#'     lopar <- min(val)
#'     hipar <- max(val)
#'     names(lopar) <- names(hipar) <- nm
#'     fitted(runlist(decrease = update(mod1, newpars = lopar),
#'                    increase = update(mod1, newpars = hipar)))
#'   }, SIMPLIFY = FALSE)
#' 
#' xyplot.list(parsims, superpose = TRUE, layout = c(1,NA),
#'             main = "Simple parameter perturbation example") +
#'   layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
#' 
    awbm = awbm.ranges(),


#' Simple degree day factor snow model
#' 
#' Simple degree day factor snow model coupled with IHACRES CMD soil moisture
#' model.
#' 
#' SWE snow water equivalent
#' 
#' ISWE water equivalent of ice in the snowpack
#' 
#' LSWE liquid water retained in the snowpack
#' 
#' @aliases snow snow.sim
#' @param DATA a \code{\link{ts}}-like object with named columns: \describe{
#' \item{list("P")}{ time series of areal rainfall depths, usually in mm. }
#' \item{list("E")}{ time series of potential evapo-transpiration, or more
#' typically, temperature as an indicator of this. } }
#' @param Tmax temperature threshold for rain, all rain is liquid above this
#' threshold.
#' @param Tmin temperature threshold for rain, all rain is snow below this
#' threshold.
#' @param Tmelt temperature threshold for snowmelt and freezing in the
#' snowpack.
#' @param kd degree day factor for snowmelt.
#' @param kf degree day factor for freezing.
#' @param rcap retention parameter for liquid water capacity of snowpack.
#' @param cr correction factor for rainfall.
#' @param cs correction factor for snowfall.
#' @param LSWE_0,ISWE_0 initial values of state variables.
#' @param \dots parameters for the \link{IHACRES.CMD.model}.
#' @param return_state to return state variables as well as the effective
#' rainfall.
#' @return \code{snow.sim} returns the modelled time series of effective
#' rainfall, or if \code{return_state = TRUE}, a multi-variate time series with
#' named columns \code{U} (effective rainfall), \code{SWE} (snow water
#' equivalent) and \code{TF}, as well as the CMD state variables.
#' @author Coded in R by Jarkko Koskela @@tkk.fi 2010-02-26.
#' 
#' Converted to C by Felix Andrews \email{felix@@nfrac.org}.
#' @seealso \code{\link{hydromad}(sma = "snow")} to work with models as objects
#' (recommended).
#' @references Kokkonen T., Jakeman A.J, Koivusalo.H, Norton.J.: COMPUTATIONAL
#' METHODS FOR WATER RESOURCE ASSESSMENTS: AN EXERCISE KIT Educational Series
#' on Modelling and Software iEMSs International Modelling and Software Society
#' Available through www.iemss.org
#' @keywords models
#' @examples
#' 
#' ## view default parameter ranges:
#' str(hydromad.options("snow"))
#' 
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "snow", routing = "expuh")
#' mod0
#' 
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, Tmax = 15, Tmin = 5, cr = 1, cs = 1, 
#'                kd = 3, kf = 1, rcap = 0.5,
#'                d = 200, f = 0.5, e = 0.1, tau_s = 10)
#' 
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[,1:2], snow = testQ))
#' 
#' ## show effect of increase/decrease in each parameter
#' parlist <- list(Tmax = c(10, 20), Tmin = c(0, 10),
#'                 cr = c(0.5, 2), cs = c(0.5, 2),
#'                 kd = c(2, 5), kf = c(0, 2), rcap = c(0, 1))
#' parsims <- mapply(val = parlist, nm = names(parlist),
#'   FUN = function(val, nm) {
#'     lopar <- min(val)
#'     hipar <- max(val)
#'     names(lopar) <- names(hipar) <- nm
#'     fitted(runlist(decrease = update(mod1, newpars = lopar),
#'                    increase = update(mod1, newpars = hipar)))
#'   }, SIMPLIFY = FALSE)
#' 
#' xyplot.list(parsims, superpose = TRUE, layout = c(1,NA),
#'             strip = FALSE, strip.left = TRUE,
#'             main = "Simple parameter perturbation example") +
#'   latticeExtra::layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
#' 
    snow = snow.ranges(),


#' Simple time-varying runoff proportion
#' 
#' Simple time-varying runoff proportion. Rainfall is scaled by the runoff
#' coefficient estimated in a moving window.  This SMA uses streamflow data, so
#' can not be used for prediction.
#' 
#' 
#' @aliases runoffratio runoffratio.sim absorbScale.hydromad.runoffratio
#' @param DATA time-series-like object with columns \code{P} (precipitation)
#' and \code{Q} (streamflow).
#' @param width width of the time window (in time steps) in which to estimate
#' the runoff coefficient.
#' @param kernel type of window used to estimate the runoff coefficient: 1 is
#' rectangular, 2 is triangular-weighted, 3 is Gaussian-like.
#' @param sides 2 for time-centered estimates, 1 for estimates using data
#' backward in time only.
#' @param rrthresh a theshold value of the runoff ratio, below which there is
#' no effective rainfall.
#' @param qlag number of time steps to lag the streamflow (relative to
#' rainfall) before estimating the runoff coefficient.
#' @param scale constant multiplier of the result, for mass balance.  If this
#' parameter is set to \code{NA} (as it is by default) in
#' \code{\link{hydromad}} it will be set by mass balance calculation.
#' @param return_state ignored.
#' @return the simulated effective rainfall, a time series of the same length
#' as the input series.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "runoffratio")} to work with models as
#' objects (recommended).
#' @keywords models
#' @examples
#' 
#' ## view default parameter ranges:
#' str(hydromad.options("runoffratio"))
#' 
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "runoffratio", routing = "expuh")
#' mod0
#' 
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, width = 30, rrthresh = 0.2, tau_s = 10)
#' 
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[,1:2], runoffratio = testQ))
#' 
#' ## show effect of increase/decrease in each parameter
#' parRanges <- list(width = c(10, 180), qlag = c(-30, 30),
#'                   rrthresh = c(0,0.5))
#' parsims <- mapply(val = parRanges, nm = names(parRanges),
#'   FUN = function(val, nm) {
#'     lopar <- min(val)
#'     hipar <- max(val)
#'     names(lopar) <- names(hipar) <- nm
#'     fitted(runlist(decrease = update(mod1, newpars = lopar),
#'                    increase = update(mod1, newpars = hipar)))
#'   }, SIMPLIFY = FALSE)
#' 
#' xyplot.list(parsims, superpose = TRUE, layout = c(1,NA),
#'             main = "Simple parameter perturbation example") +
#'   latticeExtra::layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
#' 
    runoffratio = runoffratio.ranges(),


#' Typical initial model used in Data-Based Mechanistic modelling.
#' 
#' Typical initial model used in Data-Based Mechanistic modelling.  Rainfall is
#' scaled by corresponding streamflow values raised to a power.  This SMA uses
#' streamflow data, so can not be used for prediction.
#' 
#' 
#' @aliases dbm dbm.sim absorbScale.hydromad.dbm
#' @param DATA time-series-like object with columns \code{P} (precipitation)
#' and \code{Q} (streamflow).
#' @param power power to apply to streamflow values.
#' @param qlag number of time steps to lag the streamflow (relative to
#' rainfall) before multiplication.
#' @param scale constant multiplier of the result, for mass balance.  If this
#' parameter is set to \code{NA} (as it is by default) in
#' \code{\link{hydromad}} it will be set by mass balance calculation.
#' @param return_state ignored.
#' @return the simulated effective rainfall, a time series of the same length
#' as the input series.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "dbm")} to work with models as objects
#' (recommended).
#' @keywords models
#' @examples
#' 
#' ## view default parameter ranges:
#' str(hydromad.options("dbm"))
#' 
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "dbm", routing = "expuh")
#' mod0
#' 
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, power = 0.5, qlag = 0, tau_s = 10)
#' 
#' xyplot(cbind(HydroTestData, dbm.Q = predict(mod1)))
#' 
#' ## show effect of increase/decrease in each parameter
#' parRanges <- list(power = c(0.01, 0.9), qlag = c(-1, 2))
#' parsims <- mapply(val = parRanges, nm = names(parRanges),
#'   FUN = function(val, nm) {
#'     lopar <- min(val)
#'     hipar <- max(val)
#'     names(lopar) <- names(hipar) <- nm
#'     fitted(runlist(decrease = update(mod1, newpars = lopar),
#'                    increase = update(mod1, newpars = hipar)))
#'   }, SIMPLIFY = FALSE)
#' 
#' xyplot.list(parsims, superpose = TRUE, layout = c(1,NA),
#'             main = "Simple parameter perturbation example") +
#'   layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
#' 
    dbm = dbm.ranges(),


#' Power law transfer function models
#' 
#' A power-law form of unit hydrograph (transfer function).
#' 
#' 
#' The power law form of the unit hydrograph is:
#' 
#' \deqn{H = 1 / (1 + (t/a)^{b/c}) ^ c}
#' 
#' where H is the fraction of peak flow, t is the time since peak, and a, b and
#' c are parameters.
#' 
#' From Croke (2006):
#' 
#' Parameter a is the value of t (time since peak) at which the ordinate of the
#' asymptote \eqn{(t/a)^(-b)} has a value of 1, b determines the persistence of
#' the flow response and c defines the shape of the response curve near its
#' peak. The c parameter appears twice in order to reduce interaction between
#' the b and c parameters (in this form, the c parameter only influences the
#' curvature near t = a, and doesn't influence the asymptote, which is
#' determined solely by the b parameter). The time for H to decrease to 0.5 is
#' \eqn{a(2^(1/c) - 1)^(c/b)}. While this is a three parameter model, for
#' \eqn{t >> a} only the b parameter is significant. Since the value of the a
#' parameter is typically significantly less than one (see Table 1) the
#' recession curve can be written as
#' 
#' \deqn{H = (t_r / t)^b}
#' 
#' where \eqn{t_r} is some reference time (\eqn{t_r >> a}) at which the
#' hydrograph profile has been normalized. Thus the remaining two parameters (a
#' and c) only influence the response curve near the event peak, and [the
#' equation above] can be taken as a single parameter recession model.
#' 
#' @aliases powuh powuh.sim ssg.powuh normalise.powuh
#' @param U input time series.
#' @param delay lag (dead time) between input and response, in time steps.
#' @param a the time for flow to drop by half after a peak, if \code{c = 1}.
#' See Details.
#' @param b persistence of the flow response; defines the recession curve tail.
#' @param c curvature at half-peak point.
#' @param init initial flow value(s) used in convolution filter.
#' @param uhsteps number of time steps to use in approximating the unit
#' hydrograph convolution filter.
#' @param na.action function to remove missing values, e.g.
#' \code{\link[=na.omit.ts]{na.omit}}.
#' @param epsilon values smaller than this will be set to zero.
#' @return the model output as a \code{\link{ts}} object, with the same
#' dimensions and time window as the input \code{U}.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{expuh}} \code{\link{armax}}
#' @references Croke, B.F.W. (2006). A technique for deriving an average event
#' unit hydrograph from streamflow-only data for ephemeral quick-flow-dominant
#' catchments. \emph{Advances in Water Resources} 29, pp. 493--502.
#' @keywords ts
#' @examples
#' 
#' U  <- ts(c(1, rep(0, 99)))
#' xyplot(cbind("a = 5" = powuh.sim(U, a = 5),
#'              "& b = 2" = powuh.sim(U, a = 5, b = 2),
#'              "& c = 2" = powuh.sim(U, a = 5, c = 2)),
#'        superpose = TRUE)
#' 
    powuh = powuh.ranges(),
    order = c(n = 1, m = 0),
    delay = NA,
    max.delay = 10,
    warmup = 100,
    normalise = TRUE,
    loglik = function(Q, X, ...) -0.5 * sum((Q - X)^2, na.rm = TRUE),
    objective = ~ 0.7 * hmadstat("r.sq.sqrt")(Q, X) + 0.3 * .(hmadstat("r.sq.monthly", DATA = DATA))(Q, X),
    summary.stats = c("rel.bias", "r.squared", "r.sq.sqrt", "r.sq.log"),
    prefilter = list(0.9, 0.98, 0.5),
    sriv.iterations = 12,
    sriv.epsilon = 1e-3,
    riv.noise.order = NULL,
    inverse.fit.method = "sriv",
    inverse.iterations = 30,
    inverse.rel.tolerance = 1 / 5000,
    inverse.par.epsilon = 0.001,
    sce.control = list(),
    de.control = list(itermax = 1000 / 50),
    dream.control = list(),
    cmaes.control = list(),
    nsga2.control = list(),
    dds.control = list(
      logfile = NULL,
      projectfile = NULL,
      load_projectfile = "no"
    ),
    fit.samples = 100,
    optim.method = "PORT",
    optim.control = list(
      reltol = 1e-6, maxit = 150,
      trace = if (interactive()) 1 else 0, REPORT = 10
    ),
    nlminb.control = list(
      eval.max = 200, iter.max = 150, abs.tol = 0,
      trace = if (interactive()) 10 else 0, rel.tol = 1e-6
    ),
    sim.epsilon = 1e-5,
    quiet = FALSE,
    trace = FALSE,
    catch.errors = TRUE,
    catch.errors.optim = TRUE,
    pure.R.code = FALSE,
    parallel = list(
      update.runlist = "none",
      objFunVal.runlist = "none",
      evalPars = "none",
      evalParsTS = "none",
      crossValidate = "none",
      paretoObjectivesVaryWeights = "none"
    )
  )
}

## code below copied from lattice

hydromad.getOption <- function(name) {
  .HydromadEnv$options[[name]]
}



#' User default settings for hydromad
#' 
#' A basic user settings facility, like \code{\link{options}} and
#' \code{\link{lattice.options}}.
#' 
#' These functions are direct copies of the lattice equivalents: see
#' \code{\link{lattice.options}}.
#' 
#' The available options can be seen with \code{str(hydromad.options())}.  Many
#' of these simply provide defaults for corresponding arguments to the
#' \code{\link{hydromad}} function.
#' 
#' @aliases hydromad.options hydromad.getOption
#' @param name character giving the name of a setting.
#' @param ...  new options can be defined, or existing ones modified, using one
#' or more arguments of the form 'name = value' or by passing a list of such
#' tagged values.  Existing values can be retrieved by supplying the names (as
#' character strings) of the components as unnamed arguments.
#' @seealso \code{\link{hydromad}}
#' @keywords programming
#' @examples
#' 
#' oopt <- hydromad.options()
#' str(oopt)
#' 
#' ## reset
#' hydromad.options(oopt)
#' 
#' @export hydromad.options
hydromad.options <- function(...) {
  ## this would have been really simple if only form allowed were
  ## lattice.options("foo", "bar") and
  ## lattice.options(foo=1, bar=2). But it could also be
  ## lattice.options(foo=1, "bar"), which makes some juggling necessary

  new <- list(...)
  if (is.null(names(new)) && length(new) == 1 && is.list(new[[1]])) new <- new[[1]]
  old <- .HydromadEnv$options

  ## if no args supplied, returns full options list
  if (length(new) == 0) {
    return(old)
  }

  nm <- names(new)
  if (is.null(nm)) {
    return(old[unlist(new)])
  } ## typically getting options, not setting
  isNamed <- nm != "" ## typically all named when setting, but could have mix
  if (any(!isNamed)) nm[!isNamed] <- unlist(new[!isNamed])

  ## so now everything has non-"" names, but only the isNamed ones should be set
  ## everything should be returned, however

  retVal <- old[nm]
  names(retVal) <- nm
  nm <- nm[isNamed]

  ## this used to be

  ## modified <- updateList(retVal[nm], new[nm])
  ## .LatticeEnv$lattice.options[names(modified)] <- modified

  ## but then calling lattice.options(foo = NULL) had no effect
  ## because foo would be missing from modified.  So, we now do:

  updateList <- function(x, val) {
    if (is.null(x)) x <- list()
    utils::modifyList(x, val)
  }
  .HydromadEnv$options <- updateList(old, new[nm])

  ## return changed entries invisibly
  invisible(retVal)
}
