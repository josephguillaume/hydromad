\name{gr7j}
\encoding{UTF-8}
\alias{gr7j}
\alias{gr7j.sim}
\alias{gr4jrouting}
\alias{gr4jrouting.sim}
\title{ GR7J rainfall runoff model }
\description{
  GR7J model : modified version of GR4J to allow the sensitivity of storage dependant
  evaporation and streamflow production to differ effectively creating "catchment memory". 
}
\usage{
gr7j.sim(DATA, x1, x5, x6, x7, etmult = 1, S_0 = 0.5,
  return_state = FALSE,transformed=FALSE)

gr7jrouting.sim(U, x2, x3, x4, R_0 = 0, split = 0.9,
                return_components = FALSE,
                epsilon = hydromad.getOption("sim.epsilon"),
		transformed=FALSE)
}
\arguments{
  \item{DATA}{time-series-like object with columns P (precipitation,
    mm) and E (potential evapo-transpiration, mm).}
  \item{U}{effective rainfall series.}
  \item{x1}{ maximum capacity of the production store (mm). }
  \item{x2}{ groundwater exchange coefficient (mm). }
  \item{x3}{ one day ahead maximum capacity of the routing store (mm). }
  \item{x4}{ time base of unit hydrograph UH1 (time steps). }
  \item{x5}{ sensitivity for P_s. }
  \item{x6}{ magnitude factor for E_s. }
  \item{x7}{ shape factor for E_s. }
  \item{etmult}{Multiplier for the \code{E} input data.}
  \item{S_0}{Initial soil moisture level as fraction of \code{x1}.}
  \item{R_0}{Initial groundwater reservoir level as fraction of
\code{x3}.}
  \item{split}{ Fraction to go into quick flow routing, usually fixed at 0.9. }
  \item{return_state}{
    to return the series U, S (storage) and ET (evapotranspiration).
  }
  \item{return_components}{
    to return the series Xr, Xd and R (reservoir level).
  }
  \item{epsilon}{
    values smaller than this in the output will be set to zero.
  }
  \item{transformed}{transform parameters before use to improve identifiability. They can be untransformed using \code{\link{gr4j.transformpar}}}
}
\details{
  The default parameter ranges were taken from supporting information 
  https://doi.org/10.1002/hyp.13282.
}
\value{
  the simulated effective rainfall, a time series of the same length as
  the input series.
}
\references{
Grigg AH, Hughes JD. "Nonstationarity driven by multidecadal change 
in catchment groundwater storage: A test of modifications to a common
 rainfall–runoff model." \emph{Hydrological Process.} 2018;1–14.
\url{https://doi.org/10.1002/hyp.13282}
}
\author{ Justin Hughes \email{justin.hughes@csiro.au} and Lynn Seo \email{lynn.seo@gmail.com} 
}
\seealso{
  \code{\link{?gr8j.sim}()}
}
\examples{
dat <-
  cbind(P = c(0,0,0,0,0,0.04,0.59,0.03,0.01,0.16,0.37,8.76,2.65,
              0.05,0.02,0.02,0.38,0.00,0.02,0.46,4.46,7.71,5.71,0.79,1.33),
        E = c(0,0,0,0,0,0.24,0.24,0.24,0.24,0.24,0.25,0.25,0.26,
              0.27,0.28,0.32,0.33,0.34,0.35,0.36,0.36,0.37,0.37,0.38,0.38))
datz <- zoo(dat, as.Date("2000-01-01") + 1:nrow(dat))
modz <- hydromad(datz, sma = "gr7j", routing = "gr7jrouting",
                 x1 = 665, x2 = 1.18, x3 = 90, x4 = 3.8, 
                 x5=1, x6=1, x7= 1, S_0 = 0.6, R_0 = 0.7)
xyplot(predict(modz, return_state = TRUE, return_components = TRUE),
       strip = FALSE, strip.left = TRUE)
}
\keyword{models}

