\name{IHACRES.CWI.model}
\alias{IHACRES.CWI.model}
\alias{cwi}
\alias{cwi.sim}
\alias{absorbScale.hydromad.cwi}
\title{ IHACRES Catchment Wetness Index (CWI) model }
\description{
  The Catchment Wetness Index (CWI) effective rainfall model for
  IHACRES. This is the classic model of Jakeman and Hornberger (1993),
  with the extensions to ephemeral catchments of Ye et al. (1997).  
}
\usage{
cwi.sim(DATA,
	tw, f = 0, scale,
        l = 0, p = 1,
        t_ref = hydromad.getOption("cwi")$t_ref,
        s_0 = 0,
        return_state = FALSE)
}
\arguments{
  \item{DATA}{
    a \code{\link{ts}}-like object with named columns:
    \describe{
      \item{\code{P}}{
	time series of areal rainfall depths, usually in mm.
      }
      \item{\code{E}}{
	time series of potential evapo-transpiration, or more typically,
	temperature as an indicator of this.
	Can be omitted if \code{f = 0}.
      }
    }
  }
  \item{tw}{
    drying rate at reference temperature (\eqn{\tau_\omega}{tw}).
    This is a \emph{time constant}, the number of time steps to reduce
    to a fraction \eqn{1/e \approx 37\%}. See definition below.
  }
  \item{f}{
    temperature dependence of drying rate. See definition below.
    The case of \code{f=0} describes an invariant drying rate model,
    in which case the input data \code{E} is not required.
  }
  \item{scale}{
    mass balance term (\var{c} in the literature).
    If this parameter is set to \code{NA} (as it is by default)
    in \code{\link{hydromad}} it will be set by mass balance calculation. 
  }
  \item{l}{
    moisture threshold for producing flow (in units of \var{s}).
    This can be used together with \code{p} for ephemeral rivers.
  }
  \item{p}{
    power on soil moisture (above the threshold \code{l}).
  }
  \item{t_ref}{
    reference temperature in units of \code{E} (traditionally 20 deg. C).
    This is not a parameter; it simply transforms \code{tw}
    (scaling \code{tw} by \code{exp(0.062 * t_ref * f)}.
  }
  \item{s_0}{
    starting value for soil moisture index \var{s}.
  }
  \item{return_state}{
    to return state variables as well as the effective rainfall.
  }
}
\details{
  The IHACRES model with an antecedent precipitation index was introduced by
  Jakeman et al. (1990), based on the Bedford-Ouse model of Whitehead
  et al. (1979). This slightly more physics-based version with a
  Catchment Wetness Index (CWI) was developed by Jakeman and Hornberger
  (1993). It is a metric-type model, where rainfall effectiveness is
  proportional to a simple antecedent moisture index, and the output is
  scaled to enforce mass balance.
  
	The effective rainfall at each time step is proportional to
	rainfall, scaled by a soil moisture index \var{s}:
	\deqn{U_t = c \cdot s_t \cdot P_t}{
	      U[t] = c * s[t] * P[t]}

	Or, if the parameters \code{l} and \code{p} for ephemeral rivers
	are used (after Ye et al., 1997):
	\deqn{U_t = (c (s_t - l))^p \cdot P_t}{
	      U[t] = (c * (s[t] - l))^p * P[t]}

	The soil moisture index \var{s} is calculated by a filter applied to the rainfall,
	where the drying rate is defined by a \emph{time constant} \eqn{\tau_{\omega,~t}}{tw[t]}:
	\deqn{s_t = (1 - 1 / \tau_{\omega,t}) s_{t-1} + P_t}{
	      s[t] = (1 - 1 / tw[t]) * s[t-1] + P[t]}

	If \code{f = 0} then the drying time constant is equal to the value of \code{tw}.
	Otherwise the drying rate varies over time according to the input data \code{E}:
	\deqn{\tau_{\omega,t} = \tau_\omega \exp(-0.062 f E_t)}{
	      tw[t] = tw * \exp(- 0.062 * f * E[t])}

	Note that the drying rate and effective rainfall are bounded below by 0,
	a step omitted in the equations above.
}
\value{
  \code{cwi.sim} returns the modelled time series of effective rainfall,
  or if \code{return_state = TRUE}, a multi-variate time series with named
  columns \code{U} (effective rainfall), \code{s} (index of soil
  moisture, \var{s}) and \code{w} (the recession rate of \var{s}, i.e.
  \eqn{(1 - 1 / \tau_{\omega,t})}{
    (1 - 1 / tw[t])}.
}
\references{
  Jakeman, A. J., and G. M. Hornberger (1993),
  How much complexity is warranted in a rainfall-runoff model?,
  \emph{Water Resources Research}, 29: 2637-2649.

  Jakeman, A.J., I.G. Littlewood, and P.G. Whitehead (1990),
  Computation of the instantaneous unit hydrograph and identifiable component flows
  with application to two small upland catchments,
  \emph{Journal of Hydrology}, 117: 275-300.

  Ye, W., B.C. Bates, N.R. Viney, M. Sivapalan and A.J. Jakeman (1997),
  Performance of conceptual rainfall-runoff models in low-yielding ephemeral catchments,
  \emph{Water Resources Research}, 33: 153-16.
}
\author{ Felix Andrews \email{felix@nfrac.org} }
\seealso{
  \code{\link{hydromad}(sma = "cwi")} to work with models as objects (recommended).
}
\examples{
## view default parameter ranges:
str(hydromad.options("cwi"))

data(Canning)
x <- cwi.sim(Canning[1:1000,], tw = 162, f = 2, l = 300,
             t_ref = 0, scale = 0.000284, return_state = TRUE)
xyplot(x)

data(HydroTestData)
mod0 <- hydromad(HydroTestData, sma = "cwi", routing = "expuh")
mod0

## simulate with some arbitrary parameter values
mod1 <- update(mod0, tw = 32, f = 2, scale = 0.01, tau_s = 10)

## plot results with state variables
testQ <- predict(mod1, return_state = TRUE)
xyplot(cbind(HydroTestData[,1:2], cwi = testQ))

## show effect of increase/decrease in each parameter
parRanges <- list(tw = c(0, 100), f = c(0, 8))
parsims <- mapply(val = parRanges, nm = names(parRanges),
  FUN = function(val, nm) {
    lopar <- min(val)
    hipar <- max(val)
    names(lopar) <- names(hipar) <- nm
    fitted(runlist(decrease = update(mod1, newpars = lopar),
                   increase = update(mod1, newpars = hipar)))
  }, SIMPLIFY = FALSE)

xyplot.list(parsims, superpose = TRUE, layout = c(1,NA),
            main = "Simple parameter perturbation example") +
  latticeExtra::layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
}
\keyword{ models }
