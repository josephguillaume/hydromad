#' ARMAX Transfer Function models
#' 
#' @name armax
#' 
#' @description ARMAX linear transfer functions with a single input and single output
#' series. Can be used as a general Unit Hydrograph transfer function, defined
#' by Auto-Regressive and Moving Average coefficients.
#' 
#' The transfer function used here, with input \var{u} and output \var{x} is:
#' \deqn{x[t] = a_1 x[t-1] + \ldots + a_n x[t-n] + }{ x[t] = a[1] x[t-1] + ...
#' + a[n] x[t-n] + b[0] u[t-d] + ... + b[m] u[t-m-d]}\deqn{ b_0 u[t-\delta] +
#' \ldots + b_m u[t-m-\delta]}{ x[t] = a[1] x[t-1] + ... + a[n] x[t-n] + b[0]
#' u[t-d] + ... + b[m] u[t-m-d]}
#' 
#' and the \emph{order} is denoted \eqn{(n, m)}, with delay \eqn{\delta}{d}.
#' 
#' @aliases armax armax.sim ssg.armax normalise.armax
#' @param U input time series.
#' @param a_1,a_2,a_3,b_0,b_1,b_2,b_3 ARMAX coefficients. Auto-regressive terms
#' begin with \code{a} and moving average terms begin with \code{b}. See
#' Details section.
#' @param pars the ARMAX coefficients as a named vector. If this is given, it
#' will over-ride the named parmameter arguments. Any number of terms can be
#' given here, it is not limited to the named arguments.
#' @param delay lag (dead time) between input and response, in time steps.
#' @param init initial values for the autoregressive filter.
#' @param na.action function to remove missing values, e.g.
#' \code{\link[=na.omit.ts]{na.omit}}.
#' @param epsilon values smaller than this will be set to zero.
#' @param return_components whether to return exponential component time
#' series.  If \code{TRUE}, the parameters will be converted to an exponential
#' components formulation, and passed to \code{\link{expuh.sim}}. This may fail
#' in some cases.
#' @param theta the parameters as a named vector.
#' @return the model output as a \code{\link{ts}} object, with the same
#' dimensions and time window as the input \code{U}.  If
#' \code{return_components = TRUE}, it will have multiple columns named
#' \code{Xs}, \code{Xq} and, if relevant, \code{X3}.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{armax.sriv.fit}}, \code{\link{arima}}
#' @references Jakeman, A.J., I.G. Littlewood, and P.G. Whitehead (1990),
#' Computation of the instantaneous unit hydrograph and identifiable component
#' flows with application to two small upland catchments, \emph{Journal of
#' Hydrology}, 117: 275-300.
#' @keywords ts
#' @examples
#' 
#' data(HydroTestData)
#' fit <- hydromad(HydroTestData, routing = "armax",
#'                 rfit = list("ls", order = c(n = 2, m = 1)))
#' pars <- coef(fit)
#' pars
#' 
#' xyplot(armax.sim(HydroTestData[,"P"], pars = pars))
#' 
NULL





#' Rainfall and streamflow for Bingham River Trib at Ernies Catchment.
#' 
#' Daily rainfall and streamflow for Bingham River Trib at Ernies Catchment
#' (Western Australia), from 1974-05-18 to 2008-11-02.  The catchment area is
#' 2.68 square kilometers.
#' 
#' Salmon Brook is located in the low-relief Darling Range of southwestern
#' Western Australia.
#' 
#' \describe{ \item{Rainfall (P)}{ Daily rainfall (mm/day).  Rain gauge station
#' ID 509249 "BINGHAM RIVER TRIB @ ERNIES CATCHMENT".  Latitude -33.2921;
#' Longitude 116.4451.  } \item{Streamflow (Q)}{ Daily mean streamflow
#' (mm/day).  Stream gauge ID 612008 "BINGHAM RIVER TRIB @ ERNIES CATCHMENT".
#' Latitude -33.2939; Longitude 116.4449.  } \item{Temperature (E)}{ Mean
#' Maximum Temperature Climate Data.  Product code: IDCJAC0002 reference:
#' 00232595.
#' 
#' Bureau of Meteorology station number: 9534 Station name: DONNYBROOK Latitude
#' -33.57; Longitude 115.82; Altitude 63m.  } }
#' 
#' @name BinghamTrib
#' @docType data
#' @format A \code{\link{zoo}} object, of class \code{c("zooreg", "zoo")}.  It
#' is a regular time series indexed by days, in \code{Date} format.
#' 
#' There are three columns, \code{P} (rainfall, mm/day) and \code{Q}
#' (streamflow, mm/day).  \code{E} (temperature, degrees C).
#' @seealso \code{\link{SalmonBrook}}
#' @source Water INformation (WIN) database - discrete sample data.
#' [2009-07-09]. Department of Water, Water Information Provision section,
#' Perth Western Australia.
#' 
#' Hydstra database - time-series data.  [2009-07-09].  Department of Water,
#' Water Information Provision section, Perth Western Australia.
#' 
#' Copyright (c) State of Western Australia (Department of Water).
#' 
#' Temperature: Copyright (c) Commonwealth of Australia. Created on Tue 07 Jul
#' 2009 05:48:38 AM GMT from Climate Data Online, Bureau of Meteorology.
#' http://www.bom.gov.au/climate/averages
#' @keywords datasets
#' @examples
#' 
#' data(BinghamTrib)
#' summary(BinghamTrib)
#' xyplot(BinghamTrib)
#' 
NULL





#' Rainfall, streamflow and potential evaporation data for Canning River at
#' Scenic Drive.
#' 
#' Daily rainfall, streamflow and potential evaporation for Canning River at
#' Scenic Drive upstream of Canning Dam (Western Australia), from 1977-01-01 to
#' 1987-12-31.  The catchment area is 517 square kilometers.
#' 
#' The Canning River is a tributary of the Swan River, and is located South
#' East of Perth, Western Australia.
#' 
#' \describe{ \item{Rainfall (P)}{ Daily areal rainfall (mm/day).  Origin
#' unknown, probably spatial interpolation by Barry Croke.  } \item{Streamflow
#' (Q)}{ Daily mean streamflow (mm/day).  Stream gauge ID 616024 "Canning River
#' at Scenic Drive".  Latitude -33.4176; Longitude 115.9817.  } \item{Potential
#' Evapotranspiration (E)}{ Origin Unknown.  } }
#' 
#' @name Canning
#' @docType data
#' @format A \code{\link{zoo}} object, of class \code{c("zooreg", "zoo")}.  It
#' is a regular time series indexed by days, in \code{Date} format.
#' 
#' There are three columns, \code{P} (areal rainfall, mm/day) and \code{Q}
#' (streamflow, mm/day).  \code{E} (potential evapotranspiration, mm/day).
#' @source Talk to Barry Croke...
#' @keywords datasets
#' @examples
#' 
#' data(Canning)
#' summary(Canning)
#' xyplot(Canning)
#' 
NULL





#' Unstable/unoptimised version of IHACRES Catchment Moisture Deficit (CMD)
#' model
#' 
#' Unoptimised version of CMD model. \code{g} is directly specified, and
#' therefore highly correlated with \code{d}.  Other anything than for
#' demonstration purposes, \code{\link{cmd}} should be used instead.
#' 
#' See \code{\link{cmd.sim}} for details.
#' 
#' This version is modified so that \code{g} is specified directly instead of
#' being calculated as \code{g=f*d}
#' 
#' @aliases cmd_unstable cmd_unstable.sim
#' @param DATA a \code{\link{ts}}-like object with named columns: \describe{
#' \item{list("P")}{ time series of areal rainfall depths, usually in mm. }
#' \item{list("E")}{ time series of potential evapo-transpiration, or more
#' typically, temperature as an indicator of this. } }
#' @param g CMD stress threshold
#' @param e temperature to PET conversion factor.
#' @param d CMD threshold for producing flow.
#' @param shape defines form of the \eqn{dU/dP} relationship: \code{shape = 0}
#' is the linear form, \code{shape = 1} is the trigonometric form, and
#' \code{shape > 1} is the power form.
#' @param M_0 starting CMD value.
#' @param return_state to return state variables as well as the effective
#' rainfall.
#' @return \code{cmd_unstable.sim} returns the modelled time series of
#' effective rainfall, or if \code{return_state = TRUE}, a multi-variate time
#' series with named columns \code{U} (effective rainfall), \code{CMD} and
#' \code{ET} (evapo-transpiration \eqn{E_T}).
#' @note Normally compiled C code is used for simulation, but if
#' \code{return_state = TRUE} a slower implementation in R is used.
#' @author Joseph Guillaume
#' @seealso \code{\link{cmd.sim}} for preferred version.
#' @keywords models
NULL





#' Rainfall and streamflow for Cotter River at Gingera.
#' 
#' Daily rainfall and streamflow for Cotter River at Gingera (Australian
#' Capital Territory), from 1966-05-01 to 2003-06-07.  The catchment area is
#' 148 square kilometers.
#' 
#' \describe{ \item{Rainfall (P)}{ Daily areal rainfall (mm/day).
#' 
#' Derived from rain gauges operated by Bureau of Meteorology and EcoWise. An
#' area-weighted average was used, with weights determined from a long-term
#' spline-interpolated rainfall surface.  } \item{Streamflow (Q)}{ Daily mean
#' streamflow (mm/day).  Stream gauge ID 410730 "Cotter @ Gingera".  Latitude
#' -33.592; Longitude 148.822.  } \item{Temperature (E)}{ Daily Maximum
#' Temperature.  Product code: IDCJAC0002 reference: 00232595.
#' 
#' Bureau of Meteorology station number: 070014 Station name: CANBERRA AIRPORT
#' Latitude -35.305; Longitude 149.201; Altitude 578.4 m.  } }
#' 
#' @name Cotter
#' @docType data
#' @format A \code{\link{zoo}} object, of class \code{c("zooreg", "zoo")}.  It
#' is a regular time series indexed by days, in \code{Date} format.
#' 
#' There are three columns, \code{P} (areal rainfall, mm/day) and \code{Q}
#' (streamflow, mm/day).  \code{E} (temperature, degrees C).
#' @source EcoWise?
#' 
#' Temperature: Copyright (c) Commonwealth of Australia.
#' @keywords datasets
#' @examples
#' 
#' data(Cotter)
#' summary(Cotter)
#' xyplot(Cotter)
#' 
NULL





#' Exponential components transfer function models
#' 
#' A unit hydrograph (linear transfer function) defined as a system of
#' exponentially receding components. Each component is defined by its time
#' constant and fractional volume, and if there are multiple (up to 3) such
#' components they may be in a parallel and/or series configuration.
#' 
#' 
#' The \code{expuh} model is a transfer function translating an input time
#' series \var{U} into an output series \var{X}.  It describes a configuration
#' of exponentially decaying components, each defined by a recession rate
#' \eqn{\alpha} and peak response \eqn{\beta}. However, in hydrology these
#' parameters are more easily interpreted in terms of time constants \eqn{\tau}
#' (number of time steps to reduce to a fraction \eqn{1/e}, 37\%) and
#' fractional volumes \var{v}. These are directly related as:
#' 
#' \deqn{\tau = -1 / \log(\alpha)}
#' 
#' \deqn{v = \beta / (1 - \alpha)}
#' 
#' If there are two components in parallel, these are conventionally called
#' slow (\var{s}) and quick (\var{q}) flow components. The total simulated flow
#' \var{X} is the sum of these; \eqn{X[t] = X_s[t] + X_q[t]}, and:
#' 
#' \deqn{X_s[t] = \alpha_s X_s[t-1] + \beta_s U[t]} \deqn{X_q[t] = \alpha_q
#' X_q[t-1] + \beta_q U[t]}
#' 
#' Two components might also be arranged in series rather than parallel, in
#' which case:
#' 
#' \deqn{X_s[t] = \alpha_s X_s[t-1] + \beta_s U[t]} \deqn{X[t] = \alpha_q
#' X[t-1] + \beta_q X_s[t]}
#' 
#' This configuration is specified by the argument \code{series = 1}. The
#' default \code{series = 0} specifies all components to be in parallel.
#' 
#' In the case of three components, with corresponding time constants
#' \eqn{\tau_s}, \eqn{\tau_q} and \eqn{tau_3} (\code{tau_s, tau_q, tau_3}),
#' there are four possible types of configuration:
#' 
#' \describe{ \item{list("series = 0")}{ all 3 components in parallel, i.e.
#' independent flows: \var{X = s + q + 3}. In this case \code{v_q} defaults to
#' \code{1 - v_s - v_3} in order to ensure that the total volume is 1.  }
#' \item{list("series = 1")}{ one component in parallel with two in series: the
#' \code{q} component is in series with the \code{3} component, and the
#' \code{s} component is in parallel: \var{X = (q * 3) + s}. In this case
#' \code{v_q} defaults to 1.  } \item{list("series = 2")}{ two components in
#' parallel with one in series: the \code{s} and \code{q} components are in
#' parallel and the \code{3} component is in series: \var{X = 3 * (s + q)}. In
#' this case \code{v_q} defaults to \code{1 - v_s} in order to ensure that the
#' total volume of the parallel component is 1. The total volume will be 1 if
#' \code{v_3} is also 1.  } \item{list("series = 3")}{ all 3 components in
#' series: \var{X = s * q * 3}.  In this case \code{v_q} defaults to 1. The
#' total volume will be 1 if \code{v_s} and \code{v_3} are also 1.  } }
#' 
#' @aliases expuh expuh.sim ssg.expuh normalise.expuh
#' @param U input time series.
#' @param delay lag (dead time) between input and response, in time steps.
#' @param tau_s,tau_q,tau_3 time constants (\eqn{\tau}) for the exponential
#' components.
#' @param v_s,v_q,v_3 fractional volumes (\var{v}) for the exponential
#' components.
#' @param series defines the configuration of exponential components, as being
#' in parallel and/or series (for second or third order models). See Details.
#' @param loss a constant loss (or gain, if negative) term subtracted from the
#' \emph{slow} (\code{s}) component.
#' @param Xs_0,Xq_0,X3_0 initial values of the exponential components.
#' @param pars the parameters as a named vector. If this is given, it will
#' over-ride the named parmameter arguments.
#' @param return_components whether to return all component time series.
#' @param na.action function to remove missing values, e.g.
#' \code{\link[=na.omit.ts]{na.omit}}.
#' @param epsilon values smaller than this in the output will be set to zero.
#' @param theta the parameters as a named vector.
#' @return the model output as a \code{\link{ts}} object, with the same
#' dimensions and time window as the input \code{U}.  If
#' \code{return_components = TRUE}, it will have multiple columns named
#' \code{Xs}, \code{Xq} and, if relevant, \code{X3}.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{expuh.sriv.fit}} \code{\link{armax}}
#' @references Jakeman, A.J., I.G. Littlewood, and P.G. Whitehead (1990),
#' Computation of the instantaneous unit hydrograph and identifiable component
#' flows with application to two small upland catchments, \emph{Journal of
#' Hydrology}, 117: 275-300.
#' @keywords ts
#' @examples
#' 
#' data(HydroTestData)
#' mod1 <- hydromad(HydroTestData, routing = "expuh",
#'                  tau_s = 30, tau_q = 5, v_s = 0.5)
#' flowcomps <- predict(mod1, return_components = TRUE)
#' xyplot(cbind(`Slow component` = flowcomps[,"Xs"],
#'              `Total flow` = flowcomps[,1] + flowcomps[,2]),
#'        superpose = TRUE) +
#'   layer(panel.refline(h = 0))
#' 
#' U  <- ts(c(1, rep(0, 30)))
#' xyplot(cbind("tau_s = 10" = expuh.sim(U, tau_s = 10),
#'              "& tau_q = 1" = expuh.sim(U, tau_s = 10, tau_q = 1, v_s = 0.5),
#'              "&& v_s = 0.9" = expuh.sim(U, tau_s = 10, tau_q = 1, v_s = 0.9)),
#'        superpose = TRUE)
#' 
NULL





#' Exponential components transfer function models with layered slowflow stores
#' 
#' A unit hydrograph with a quickflow pathway and two layered slowflow pathways
#' modelling recharge to groundwater in order to allow modelling of long-term
#' disconnection of slowflow stores from streamflow.
#' 
#' 
#' The \code{expuh3s} model consists of a single quickflow pathway modelled as
#' an exponential store, and a slowflow pathway comprised of two layered
#' stores.
#' 
#' Each slowflow store is modelled as a \code{\link{leakyExpStore}}, which has
#' a loss term, has no flow when the store drops below a given level, and can
#' therefore model longer-term disconnection of a store from streamflow.
#' 
#' Adapted from Herron and Croke (2009):
#' 
#' The upper store, G1, receives rainfall inputs and discharges to the stream,
#' Qs and recharges the lower store. G1 has a lower limit of 0, where flow
#' ceases representing the fully 'drained' condition. Conceptually, the upper
#' store can be viewed as a perched water table, which develops in response to
#' rain and tends to be relatively short-lived, perhaps seasonal. Thus the time
#' constant, \code{tau_s}, for discharge from the 'soil' store will be
#' somewhere between that for quickflow, \code{tau_q} and the groundwater
#' discharge constant, \code{tau_g}.
#' 
#' G2 is recharged from G1 when \code{G1>G_1} and discharges to the stream
#' \code{Q_g} when \code{G2>0}. The sum of \code{Q_s} and \code{Q_g} represents
#' the total slowflow pathway. We assume that all extraction and natural
#' groundwater losses (\code{loss}) are from G2. The approach avoids the need
#' to specify a maximum capacity for either storage, but the introduction of a
#' recharge term, \code{R} between the stores adds a new parameter.
#' 
#' Recharge is represented by a constant rate \code{R} which ceases when
#' \code{G1<G_1}, diminishing linearly to that point when
#' \code{thres<G1<thres+loss}. Setting \code{G_1=0} (the default) ceases
#' recharge when flow ceases.
#' 
#' @aliases expuh3s expuh3s.sim
#' @param U input time series (units below assume ML/day)
#' @param delay lag (dead time) between input and response, in time steps.
#' @param v_s Fraction of effective rainfall that goes to groundwater
#' @param tau_q Recession coefficient for quickflow (days)
#' @param tau_s Recession coefficient for soil store (G_1) discharge (days)
#' @param tau_g Recession coefficient for groundwater store (G_2) discharge
#' (days)
#' @param R Maximum recharge from G_1 to G_2 (ML/day)
#' @param G_1 storage threshold to stop recharge (ML) (less than zero)
#' @param loss Groundwater loss (ML/day)
#' @param G_2 storage threshold to stop groundwater loss (ML) (less than zero)
#' @param Xs_0,Xq_0,X3_0 initial values of the exponential components.
#' @param pars the parameters as a named vector. If this is given, it will
#' over-ride the named parmameter arguments.
#' @param return_components whether to return all component time series.
#' @param na.action function to remove missing values, e.g.
#' \code{\link[=na.omit.ts]{na.omit}}.
#' @param epsilon values smaller than this in the output will be set to zero.
#' @return the model output as a \code{\link{ts}} object, with the same
#' dimensions and time window as the input \code{U}.  If
#' \code{return_components = TRUE}, it will have multiple columns named
#' \code{Xs}, \code{Xq} and \code{Xg}.
#' @author Joseph Guillaume \email{joseph.guillaume@@anu.edu.au}
#' @seealso \code{\link{expuh}},\code{\link{leakyExpStore}}
#' @references Herron, N.F. and B.F.W. Croke (2009). IHACRES-3S - A 3-store
#' formulation for modelling groundwater-surface water interactions. In
#' Anderssen, R.S., R.D. Braddock and L.T.H. Newham (eds) \emph{18th World
#' IMACS Congress and MODSIM09 International Congress on Modelling and
#' Simulation.} Modelling and Simulation Society of Australia and New Zealand
#' and International Association for Mathematics and Computers in Simulation,
#' July 2009, pp. 3081-3087. ISBN: 978-0-9758400-7-8.
#' \url{http://www.mssanz.org.au/modsim09/I1/herron.pdf}
#' @keywords ts
NULL





#' Support for parallelisation in hydromad
#' 
#' Hydromad allows model runs in some functions to be parallelised. The
#' parallelisation methods provided for each function depend on its
#' characteristics and are documented below and in each function's help page.
#' Parallelisation is by default turned off, as it requires the user to set it
#' up for their particular computer and to judge whether parallelisation is
#' worthwhile for their particular case.
#' 
#' A number of functions in hydromad involve performing a large number of model
#' runs that are to some extent independent of each other. Hydromad uses a
#' number of R packages to allow these models runs to occur simultaneously.
#' This usually involves splitting the job between separate 'worker' R
#' sessions.
#' 
#' However, parallelisation involves an overhead. Even on a single computer, it
#' takes extra time to communicate with workers to transfer instructions and
#' retrieve results. Parallelisation is therefore not always worthwhile, and
#' the best method for parallelising a particular analysis depends on its
#' characteristics. Do not use parallelisation if you only have a single core,
#' if you will run out of memory or the analysis already runs quickly.
#' 
#' Several settings of the parallelisation can be modified by specifying a list
#' of options to the function to be parallelised, with the argument
#' \code{parallel}, e.g.
#' 
#' \code{evalPars(pars,model,parallel=list(method="foreach",
#' packages=c("hydromad","fuse"), async=TRUE))}
#' 
#' The available options are described below.
#' 
#' \subsection{method
#' 
#' The available \code{method} are specific to each function that is
#' parallelised and documented in its respective help page. A summary of the
#' most common methods is given here.
#' 
#' \code{method="clusterApply"} uses either the built-in \code{parallel}
#' package or the \code{snow} package. Setting up the cluster requires code
#' like:
#' 
#' \code{library(parallel)}\cr \code{cl <- makeCluster(2, type="SOCK")}
#' 
#' Functions look for the object 'cl' in the global environment. The function
#' can be run on multiple cores on non-Windows machines by creating cl using
#' \code{\link{makeForkCluster}}.
#' 
#' \code{method="foreach"} allows a number of backends, e.g. the packages
#' \pkg{doParallel}, \pkg{doRedis}. Each of these packages provides a
#' backend-specific registration function that needs to be called, e.g.
#' 
#' \code{library(doParallel)}\cr \code{registerDoParallel()}
#' 
#' \code{foreach} generally incurs a higher overhead than other methods, but
#' has the advantage of flexibility. }
#' 
#' \subsection{export If the model being run depends on functions or variables
#' not defined in the hydromad package, they may need to be explicitly
#' exported. The names of the variables to be exported are specified as a
#' character vector, e.g. \code{export=c("mySMA.sim","myNewObjectiveFunction")}
#' 
#' Examples: \itemize{ \item if you have created your own soil moisture
#' accounting function \code{mySMA.sim} \item if your objective function or a
#' model you have made yourself depends on other functions or variables in the
#' global environment, e.g.
#' \code{objective=function(Q,X){myNewObjectiveFunction(Q-X)+hmadstat("r.squared")(Q,X)}}
#' } }
#' 
#' \subsection{packages If you are using models that defined in another
#' package, e.g. the \code{fuse} set of models, then the workers can be told to
#' load the package as well as hydromad by specifying, e.g.
#' \code{packages=c("hydromad","fuse")}. Make sure the package is installed on
#' all the workers. }
#' 
#' \subsection{async If the parallelisation method supports it, some functions
#' allow runs to occur in the background. Instead the function returns
#' immediately. Results can be retrieved later using methods specific to the
#' parallelisation method. Where available, this feature is enabled with
#' \code{async=TRUE} }
#' 
#' @aliases hydromad_parallelisation parallel
#' @author Joseph Guillaume
#' @seealso Functions with parallelisation: \code{\link{evalPars}},
#' \code{\link{crossValidate}}, \code{update.runlist},
#' \code{\link{objFunVal.runlist}}, \code{\link{paretoObjectivesVaryWeights}}
#' @references
#' @keywords utilities
#' @examples
#' 
#' 
NULL





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
#' ################################################################################
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
NULL





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
#' @examples
#' 
#' 
#' 
NULL





#' A simple simulated dataset for use in testing hydrological models.
#' 
#' A simple simulated dataset for use in testing hydrological models.
#' 
#' \describe{ \item{Rainfall (P)}{ a regular series of impulses (every 20 time
#' steps). Each of these pulses have a value of 5, except one which has a value
#' of 20.  } \item{Temperature (E)}{ a sine wave ranging from 0 to 30.  }
#' \item{Streamflow (Q)}{ proportional to the square of rainfall and inversely
#' to temperature, then filtered with a second-order autoregressive
#' \code{\link{filter}}.  } }
#' 
#' @format A \code{\link{zooreg}} object with 730 (365 * 2) time steps.
#' 
#' There are three columns, \code{P} (simulated areal rainfall, mm/day),
#' \code{E} (simulated temperature, degrees Celcius) and \code{Q} (simulated
#' streamflow, mm/day).
#' @keywords datasets
#' @examples
#' 
#' data(HydroTestData)
#' summary(HydroTestData)
#' xyplot(HydroTestData)
#' 
NULL





#' IHACRES Catchment Moisture Deficit (CMD) model
#' 
#' The Catchment Moisture Deficit (CMD) effective rainfall model for IHACRES.
#' It is a conceptual-type model, where input rainfall is partitioned
#' explicitly into drainage, evapo-transpiration, and changes in catchment
#' moisture.
#' 
#' The mass balance step is: \deqn{M[t] = M[t-1] - P[t] + E_T[t] + U[t]}
#' 
#' where \eqn{M} represents catchment moisture deficit (CMD), constrained below
#' by 0 (the nominal fully saturated level).  P is catchment areal rainfall,
#' \eqn{E_T} is evapo-transpiration, and U is drainage (effective rainfall).
#' All are, typically, in units of mm per time step.
#' 
#' Rainfall effectiveness (i.e. drainage proportion) is a simple
#' \emph{instantaneous} function of the CMD, with a threshold at \eqn{M = d}.
#' In the default linear form this is:
#' 
#' \deqn{\frac{\mathrm{d}U}{\mathrm{d}P} = 1 - \min(1, M/d)}{ dU/dP = 1 -
#' min(1, M/d)}
#' 
#' The trigonometric form is
#' 
#' \deqn{\frac{\mathrm{d}U}{\mathrm{d}P} = 1 - \min(1, \sin^2(\pi M / 2d))}{
#' dU/dP = 1 - min(1, sin^2(pi M / 2d))}
#' 
#' The power form is
#' 
#' \deqn{\frac{\mathrm{d}U}{\mathrm{d}P} = 1 - \min(1, (M/d)^a)}{ dU/dP = 1 -
#' min(1, (M/d)^a)} where a = 10 ^ (shape / 50)
#' 
#' The actual drainage each time step involves the integral of these relations.
#' 
#' Evapo-transpiration is also a simple function of the CMD, with a threshold
#' at \eqn{M = f d}{M = f * d}: \deqn{E_T[t] = e E[t] \min(1,
#' \exp\left(2\left(1 - \frac{M_f}{fd}\right)\right))}{ E_T[t] = e E[t] \min(1,
#' \exp(2(1 - M_f / (fd))))}
#' 
#' Note that the evapo-transpiration calculation is based on \eqn{M_f}, which
#' is the CMD after precipitation and drainage have been accounted for.
#' 
#' @aliases IHACRES.CMD.model cmd cmd.sim
#' @param DATA a \code{\link{ts}}-like object with named columns: \describe{
#' \item{list("P")}{ time series of areal rainfall depths, usually in mm. }
#' \item{list("E")}{ time series of potential evapo-transpiration, or more
#' typically, temperature as an indicator of this. } }
#' @param f CMD stress threshold as a proportion of \code{d}.
#' @param e temperature to PET conversion factor.
#' @param d CMD threshold for producing flow.
#' @param shape defines form of the \eqn{dU/dP} relationship: \code{shape = 0}
#' is the linear form, \code{shape = 1} is the trigonometric form, and
#' \code{shape > 1} is the power form.
#' @param M_0 starting CMD value.
#' @param return_state to return state variables as well as the effective
#' rainfall.
#' @return \code{cmd.sim} returns the modelled time series of effective
#' rainfall, or if \code{return_state = TRUE}, a multi-variate time series with
#' named columns \code{U} (effective rainfall), \code{CMD} and \code{ET}
#' (evapo-transpiration \eqn{E_T}).
#' @note Normally compiled C code is used for simulation, but if
#' \code{return_state = TRUE} a slower implementation in R is used.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "cmd")} to work with models as objects
#' (recommended).
#' @references Croke, B.F.W. and A.J. Jakeman (2004), A Catchment Moisture
#' Deficit module for the IHACRES rainfall-runoff model, \emph{Environmental
#' Modelling and Software}, 19(1): 1-5.
#' 
#' Croke, B.F.W. and A.J. Jakeman (2005), Corrigendum to ``A Catchment Moisture
#' Deficit module for the IHACRES rainfall-runoff model'' [Environ. Model.
#' Softw. 19 (1) (2004) 1-5], \emph{Environmental Modelling and Software},
#' 20(7): 977.
#' @keywords models
#' @examples
#' 
#' ## view default parameter ranges:
#' str(hydromad.options("cmd"))
#' 
#' data(Canning)
#' x <- cmd.sim(Canning[1:1000,], d = 200, f = 0.7, e = 0.166,
#'              return_state = TRUE)
#' xyplot(x)
#' 
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "cmd", routing = "expuh")
#' mod0
#' 
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, d = 200, f = 0.5, e = 0.1, tau_s = 10)
#' 
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[,1:2], cmd = testQ))
#' 
#' ## show effect of increase/decrease in each parameter
#' parlist <- list(d = c(50, 550), f = c(0.01, 3),
#'                 e = c(0.01, 1.5))
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
#'             main = "Simple parameter perturbation example") +
#'   layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
#' 
NULL





#' IHACRES Catchment Wetness Index (CWI) model
#' 
#' The Catchment Wetness Index (CWI) effective rainfall model for IHACRES. This
#' is the classic model of Jakeman and Hornberger (1993), with the extensions
#' to ephemeral catchments of Ye et al. (1997).
#' 
#' The IHACRES model with an antecedent precipitation index was introduced by
#' Jakeman et al. (1990), based on the Bedford-Ouse model of Whitehead et al.
#' (1979). This slightly more physics-based version with a Catchment Wetness
#' Index (CWI) was developed by Jakeman and Hornberger (1993). It is a
#' metric-type model, where rainfall effectiveness is proportional to a simple
#' antecedent moisture index, and the output is scaled to enforce mass balance.
#' 
#' The effective rainfall at each time step is proportional to rainfall, scaled
#' by a soil moisture index \var{s}: \deqn{U_t = c \cdot s_t \cdot P_t}{ U[t] =
#' c * s[t] * P[t]}
#' 
#' Or, if the parameters \code{l} and \code{p} for ephemeral rivers are used
#' (after Ye et al., 1997): \deqn{U_t = (c (s_t - l))^p \cdot P_t}{ U[t] = (c *
#' (s[t] - l))^p * P[t]}
#' 
#' The soil moisture index \var{s} is calculated by a filter applied to the
#' rainfall, where the drying rate is defined by a \emph{time constant}
#' \eqn{\tau_{\omega,~t}}{tw[t]}: \deqn{s_t = (1 - 1 / \tau_{\omega,t}) s_{t-1}
#' + P_t}{ s[t] = (1 - 1 / tw[t]) * s[t-1] + P[t]}
#' 
#' If \code{f = 0} then the drying time constant is equal to the value of
#' \code{tw}.  Otherwise the drying rate varies over time according to the
#' input data \code{E}: \deqn{\tau_{\omega,t} = \tau_\omega \exp(-0.062 f
#' E_t)}{ tw[t] = tw * \exp(- 0.062 * f * E[t])}
#' 
#' Note that the drying rate and effective rainfall are bounded below by 0, a
#' step omitted in the equations above.
#' 
#' @aliases IHACRES.CWI.model cwi cwi.sim absorbScale.hydromad.cwi
#' @param DATA a \code{\link{ts}}-like object with named columns: \describe{
#' \item{list("P")}{ time series of areal rainfall depths, usually in mm.  }
#' \item{list("E")}{ time series of potential evapo-transpiration, or more
#' typically, temperature as an indicator of this.  Can be omitted if \code{f =
#' 0}.  } }
#' @param tw drying rate at reference temperature (\eqn{\tau_\omega}{tw}).
#' This is a \emph{time constant}, the number of time steps to reduce to a
#' fraction \eqn{1/e \approx 37\%}. See definition below.
#' @param f temperature dependence of drying rate. See definition below.  The
#' case of \code{f=0} describes an invariant drying rate model, in which case
#' the input data \code{E} is not required.
#' @param scale mass balance term (\var{c} in the literature).  If this
#' parameter is set to \code{NA} (as it is by default) in
#' \code{\link{hydromad}} it will be set by mass balance calculation.
#' @param l moisture threshold for producing flow (in units of \var{s}).  This
#' can be used together with \code{p} for ephemeral rivers.
#' @param p power on soil moisture (above the threshold \code{l}).
#' @param t_ref reference temperature in units of \code{E} (traditionally 20
#' deg. C).  This is not a parameter; it simply transforms \code{tw} (scaling
#' \code{tw} by \code{exp(0.062 * t_ref * f)}.
#' @param s_0 starting value for soil moisture index \var{s}.
#' @param return_state to return state variables as well as the effective
#' rainfall.
#' @return \code{cwi.sim} returns the modelled time series of effective
#' rainfall, or if \code{return_state = TRUE}, a multi-variate time series with
#' named columns \code{U} (effective rainfall), \code{s} (index of soil
#' moisture, \var{s}) and \code{w} (the recession rate of \var{s}, i.e.
#' \eqn{(1 - 1 / \tau_{\omega,t})}{ (1 - 1 / tw[t])}.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "cwi")} to work with models as objects
#' (recommended).
#' @references Jakeman, A. J., and G. M. Hornberger (1993), How much complexity
#' is warranted in a rainfall-runoff model?, \emph{Water Resources Research},
#' 29: 2637-2649.
#' 
#' Jakeman, A.J., I.G. Littlewood, and P.G. Whitehead (1990), Computation of
#' the instantaneous unit hydrograph and identifiable component flows with
#' application to two small upland catchments, \emph{Journal of Hydrology},
#' 117: 275-300.
#' 
#' Ye, W., B.C. Bates, N.R. Viney, M. Sivapalan and A.J. Jakeman (1997),
#' Performance of conceptual rainfall-runoff models in low-yielding ephemeral
#' catchments, \emph{Water Resources Research}, 33: 153-16.
#' @keywords models
#' @examples
#' 
#' ## view default parameter ranges:
#' str(hydromad.options("cwi"))
#' 
#' data(Canning)
#' x <- cwi.sim(Canning[1:1000,], tw = 162, f = 2, l = 300,
#'              t_ref = 0, scale = 0.000284, return_state = TRUE)
#' xyplot(x)
#' 
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "cwi", routing = "expuh")
#' mod0
#' 
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, tw = 32, f = 2, scale = 0.01, tau_s = 10)
#' 
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[,1:2], cwi = testQ))
#' 
#' ## show effect of increase/decrease in each parameter
#' parRanges <- list(tw = c(0, 100), f = c(0, 8))
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
NULL





#' Runoff as rainfall to a power
#' 
#' Runoff as rainfall to a power.  This allows an increasing fraction of runoff
#' to be generated by increasingly intense/large rainfall events (for
#' \code{power > 0}).  The fraction increases up to a full runoff level at
#' \code{maxP}.
#' 
#' 
#' @aliases intensity intensity.sim absorbScale.hydromad.intensity
#' @param DATA time-series-like object with columns \code{P} (precipitation)
#' and \code{Q} (streamflow).
#' @param power power on rainfall used to estimate effective rainfall.
#' @param maxP level of rainfall at which full runoff occurs (effective
#' rainfall == rainfall).
#' @param scale constant multiplier of the result, for mass balance.  If this
#' parameter is set to \code{NA} (as it is by default) in
#' \code{\link{hydromad}} it will be set by mass balance calculation.
#' @param return_state ignored.
#' @return the simulated effective rainfall, a time series of the same length
#' as the input series.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{hydromad}(sma = "intensity")} to work with models as
#' objects (recommended).
#' @keywords models
#' @examples
#' 
#' ## view default parameter ranges:
#' str(hydromad.options("intensity"))
#' 
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "intensity", routing = "expuh")
#' mod0
#' 
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, power = 1, maxP = 200, tau_s = 10)
#' 
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[,1:2], intensity = testQ))
#' 
#' ## show effect of increase/decrease in each parameter
#' parRanges <- list(power = c(0,2), maxP = c(100,500), scale = NA)
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
NULL





#' Smoothed exponential stores
#' 
#' Multiple exponential stores with pre- and post- smoothing on 2nd and 3rd
#' stores
#' 
#' 
#' @aliases maexpuh maexpuh.sim
#' @param U placeholder
#' @param delay placeholder
#' @param tau_s placeholder
#' @param tau_q placeholder
#' @param tau_3 placeholder
#' @param v_s placeholder
#' @param v_q placeholder
#' @param v_3 placeholder
#' @param series placeholder
#' @param loss placeholder
#' @param Xs_0 placeholder
#' @param Xq_0 placeholder
#' @param X3_0 placeholder
#' @param w_s placeholder
#' @param w_3 placeholder
#' @param pars placeholder
#' @param return_components placeholder
#' @param na.action placeholder
#' @param epsilon placeholder
#' @return
#' @author Joseph Guillaume
#' @seealso \code{\link{expuh}}
#' @keywords ts
#' @examples
#' 
#' \dontrun{
#'   ## Pre and post filtered Xs
#'   Us[]=filter(U,rep(1/w_s,w_s),sides=1)
#'   Us[1:(w_s-1)] <- cumsum(U[1:(w_s-1)])/1:(w_s-1)
#'   Xstemp[] <- filter_loss(beta_s * U, alpha_s, loss = lossVal,
#'                       init = Xs_0)
#'   Xs[]=filter(Xstemp,rep(1/w_s,w_s),sides=1)
#'   Xs[1:(w_s-1)] <- cumsum(Xstemp[1:(w_s-1)])/1:(w_s-1)
#' }
NULL





#' Rainfall and streamflow for Murrindindi River at Murrindindi above Colwells.
#' 
#' Daily rainfall and streamflow for Murrindindi River at Murrindindi above
#' Colwells (Victoria, Australia), from 1975-06-08 to 1998-06-30.  The
#' catchment area is 105.9 square kilometers.
#' 
#' \describe{ \item{Rainfall (P)}{ Daily areal rainfall (mm/day).
#' 
#' Unknown...
#' 
#' Derived from rain gauges operated by Bureau of Meteorology.  An
#' area-weighted average was used, with weights determined from a long-term
#' spline-interpolated rainfall surface.  } \item{Streamflow (Q)}{ Daily mean
#' streamflow (mm/day).  Stream gauge ID 405205A "Murrindindi @ Murrindidi
#' above Colwells".  Latitude -37.412; Longitude 145.558.  } \item{Temperature
#' (E)}{ Unknown.  } }
#' 
#' @name Murrindindi
#' @docType data
#' @format A \code{\link{zoo}} object, of class \code{c("zooreg", "zoo")}.  It
#' is a regular time series indexed by days, in \code{Date} format.
#' 
#' There are three columns, \code{P} (areal rainfall, mm/day) and \code{Q}
#' (streamflow, mm/day).  \code{E} (temperature, degrees C).
#' @source Department of Sustainability and Environment (VIC).
#' 
#' Temperature: Copyright (c) Commonwealth of Australia.
#' @keywords datasets
#' @examples
#' 
#' data(Murrindindi)
#' summary(Murrindindi)
#' xyplot(Murrindindi)
#' 
NULL





#' Rainfall and streamflow for Queanbeyan River at Tinderry.
#' 
#' Daily rainfall and streamflow for Queanbeyan River at Tinderry (Australian
#' Capital Territory), from 1966-05-01 to 2003-06-07.  The catchment area is
#' 490 (BOM) or 506 (?) square kilometers.
#' 
#' \describe{ \item{Rainfall (P)}{ Daily areal rainfall (mm/day).
#' 
#' Derived from rain gauges operated by Bureau of Meteorology and EcoWise. An
#' area-weighted average was used, with weights determined from a long-term
#' spline-interpolated rainfall surface.  } \item{Streamflow (Q)}{ Daily mean
#' streamflow (mm/day).  Stream gauge ID 410734 "Queanbeyan @ Tinderry".
#' Latitude -35.615; Longitude 149.348.  } \item{Temperature (E)}{ Daily
#' Maximum Temperature.  Product code: IDCJAC0002 reference: 00232595.
#' 
#' Bureau of Meteorology station number: 070014 Station name: CANBERRA AIRPORT
#' Latitude -35.305; Longitude 149.201; Altitude 578.4 m.  } }
#' 
#' @name Queanbeyan
#' @docType data
#' @format A \code{\link{zoo}} object, of class \code{c("zooreg", "zoo")}.  It
#' is a regular time series indexed by days, in \code{Date} format.
#' 
#' There are three columns, \code{P} (areal rainfall, mm/day) and \code{Q}
#' (streamflow, mm/day).  \code{E} (temperature, degrees C).
#' @source EcoWise?
#' 
#' Temperature: Copyright (c) Commonwealth of Australia.
#' @keywords datasets
#' @examples
#' 
#' data(Queanbeyan)
#' summary(Queanbeyan)
#' xyplot(Queanbeyan)
#' 
NULL





#' Rainfall and streamflow for Salmon Brook at Salmon Catchment.
#' 
#' Daily rainfall and streamflow for Salmon Brook at Salmon Catchment (Western
#' Australia), from 1974-04-06 to 1999-03-16.  The catchment area is 0.83
#' square kilometers.
#' 
#' Salmon Brook is located in the low-relief Darling Range of southwestern
#' Western Australia. Evapotranspiration consumes about 90 per cent of the
#' annual rainfall. The region is dominated by jarrah (Eucalyptus marginata)
#' forest. Surface soils are predominantly highly permeable sands and gravels
#' (Ye et al., 1997).
#' 
#' \describe{ \item{Rainfall (P)}{ Daily rainfall (mm/day).  Rain gauge station
#' ID 509247 "SALMON BROOK @ SALMON CATCHMENT".  Latitude -33.4157; Longitude
#' 115.9846.  } \item{Streamflow (Q)}{ Daily mean streamflow (mm/day).  Stream
#' gauge ID 612011 "SALMON BROOK @ SALMON CATCHMENT".  Latitude -33.4176;
#' Longitude 115.9817.  } \item{Temperature (E)}{ Mean Maximum Temperature
#' Climate Data.  Product code: IDCJAC0002 reference: 00232595.
#' 
#' Bureau of Meteorology station number: 9534 Station name: DONNYBROOK Latitude
#' -33.57; Longitude 115.82; Altitude 63m.  } }
#' 
#' @name SalmonBrook
#' @docType data
#' @format A \code{\link{zoo}} object, of class \code{c("zooreg", "zoo")}.  It
#' is a regular time series indexed by days, in \code{Date} format.
#' 
#' There are three columns, \code{P} (rainfall, mm/day) and \code{Q}
#' (streamflow, mm/day).  \code{E} (temperature, degrees C).
#' @seealso \code{\link{BinghamTrib}}, \code{\link{YeAl97}}
#' @references Ye, W., B. C. Bates, N. R. Viney, M. Sivapalan and A. J. Jakeman
#' (1997).  Performance of conceptual rainfall-runoff models in low-yielding
#' ephemeral catchments, \emph{Water Resources Research} 33, pp. 153--166.
#' @source Water INformation (WIN) database - discrete sample data.
#' [2009-07-09]. Department of Water, Water Information Provision section,
#' Perth Western Australia.
#' 
#' Hydstra database - time-series data.  [2009-07-09].  Department of Water,
#' Water Information Provision section, Perth Western Australia.
#' 
#' Copyright (c) State of Western Australia (Department of Water).
#' 
#' Temperature: Copyright (c) Commonwealth of Australia. Created on Tue 07 Jul
#' 2009 05:48:38 AM GMT from Climate Data Online, Bureau of Meteorology.
#' http://www.bom.gov.au/climate/averages
#' @keywords datasets
#' @examples
#' 
#' data(SalmonBrook)
#' summary(SalmonBrook)
#' xyplot(SalmonBrook)
#' 
NULL





#' Rainfall and streamflow for Wye at Cefn Brwyn.
#' 
#' Hourly rainfall and streamflow for Wye at Cefn Brwyn (Wales, UK), from
#' 1987-01-01 12:00 to 1989-01-01 11:00.
#' 
#' This dataset is described and analysed in the paper cited below. It gives
#' the following introduction:
#' 
#' \dQuote{The Wye at Cefn Brwyn (the Wye) is a 10.6 km2, predominantly open
#' moorland, headwater catchment in mid-Wales in which land use is
#' predominantly sheep farming. It is one of the wettest gauged basins in
#' England and Wales; mean annual rainfall is about 2490 mm, of which about
#' 87\% leaves the catchment as streamflow (NERC, 2003). The catchment is one
#' of the Plynlimon research basins operated by the Centre for Ecology and
#' Hydrology (CEH) (e.g. Robinson and Dupeyrat, 2005). The hydrometric data
#' used in this paper for the Wye comprise 15-min flow data and hourly
#' catchment rainfall data, from which 1-, 2-, 4-, 6-, 12- and 24-hourly
#' rainfall-streamflow data sets were prepared. Streamflow for the Wye is
#' measured at a weir. The rainfall and streamflow data, extracted from the CEH
#' Plynlimon data archive, are considered to be of excellent quality.}
#' \cite{(Littlewood and Croke, 2008, p. 687)}.
#' 
#' Note: first three months of rainfall data appear to be daily averages, not
#' hourly.
#' 
#' @name Wye
#' @docType data
#' @format A \code{\link{zoo}} object, of class \code{c("zooreg", "zoo")}.  It
#' is a regular time series indexed by hours, in \code{POSIXct} format.
#' 
#' There are two columns, \code{P} (rainfall, mm / hour) and \code{Q}
#' (streamflow, mm / hour).
#' @references Littlewood, I.G. and Croke, B.F.W. (2008). Data time-step
#' dependency of conceptual rainfall-streamflow model parameters: an empirical
#' study with implications for regionalisation. \emph{Hydrological Sciences
#' Journal}, 53(4), 685-695.
#' @source \emph{Centre for Ecology and Hydrology} (CEH), Wallingford, UK, via
#' the \emph{Top-Down modelling Working Group} (TDWG) for the \emph{Prediction
#' in Ungauged Basins} (PUB) IAHS Decade (2003-2012):
#' 
#' \url{http://tdwg.catchment.org/datasets.html}
#' 
#' Thanks to Ian Littlewood for helping to organise this dataset.
#' @keywords datasets
#' @examples
#' 
#' data(Wye)
#' 
#' xyplot(Wye)
#' 
#' ## note: first three months of rainfall data are daily averages
#' xyplot(window(Wye, start = as.POSIXct("1987-01-01"),
#'                      end = as.POSIXct("1987-09-01")))
#' 
#' ## auto- and cross-correlation
#' acf(coredata(Wye[,2:1]))
#' 
NULL





#' Model performance statistics from Ye et al. 1997
#' 
#' Performance statistics from calibrations of three model structures in three
#' catchments in two data periods, assessed using absolute mean deviation, bias
#' and efficiency.
#' 
#' 
#' @name YeAl97
#' @docType data
#' @format A data frame with 36 observations on the following 7 variables.
#' \describe{ \item{list("Catchment")}{a character vector, either Canning,
#' Salmon or Stones} \item{list("calib.period")}{a character vector, either
#' "First 5Y" or "Second 5Y"} \item{list("sim.period")}{a character vector,
#' either "First 5Y" or "Second 5Y"} \item{list("Model.str")}{a factor with
#' levels \code{GSFB} \code{IHACRES} \code{LASCAM}} \item{list("A")}{a numeric
#' vector, absolute mean deviation} \item{list("B")}{a numeric vector, bias}
#' \item{list("E")}{a numeric vector, efficiency} }
#' @references Used by \code{\link{paretoTimeAnalysis}}.
#' @source Ye, W., B.C. Bates, N.R. Viney, M. Sivapalan and A.J. Jakeman
#' (1997). Performance of conceptual rainfall-runoff models in low-yielding
#' ephemeral catchments. \emph{Water Resour. Res.} 33(1): 153-166 DOI:
#' \url{http://dx.doi.org/10.1029/96wr02840}.
#' @keywords datasets
#' @examples
#' 
#' data(YeAl97)
#' YeAl97
#' 
NULL



