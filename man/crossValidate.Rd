\name{crossValidate}
\alias{crossValidate}
\title{
  Cross-validation of hydromad model specification
}
\description{
  Using a split dataset, 
  Estimate parameters using optimisation on data subsets and evaluate
  their performance on all other
}
\usage{
crossValidate(MODEL, periods, name.Model.str = paste(MODEL$sma, MODEL$routing), 
    name.Cal.objfn = "unknown", fitBy, ...)
}
\arguments{
  \item{MODEL}{
    an object of class \code{hydromad}.
  }
  \item{periods}{
    named list of start and end dates, passed to \code{\link{splitData}}
  }
  \item{name.Model.str}{
    Name to give to this model structure to allow a combined analysis
  }
  \item{name.Cal.objfn}{
    Name to give to this model identification process (e.g. name of
    objective function and/or optimisation algorithm), to allow a
    combined analysis
  }
  \item{fitBy}{
    function to estimate parameters of \code{MODEL},
    e.g. \code{\link{fitByOptim}}, \code{\link{fitBySCE}}
  }
  \item{\dots}{
    Arguments passed to \code{fitBy}
}
}
\value{
  A runlist of n*n models evaluated in each of n periods with parameters
  estimated from each of n periods, of subclass \code{crossvalidation}
}
\author{
Joseph Guillaume  
}
\seealso{
  \code{\link{paretoTimeAnalysis}}
  }
\examples{
data(Cotter)
modx <- hydromad(Cotter, sma = "cwi", routing = "expuh",
                 tau_s = c(2,100), v_s = c(0,1))
periods=list(P1=as.Date(c("1966-05-01","1976-04-30")),
             P2=as.Date(c("1976-05-1","1986-04-30")))
## Estimate parameters using single fitByOptim run
## from single initial parameter set
runs=crossValidate(modx,periods=periods,fitBy=fitByOptim,samples=1)
summary(runs)
## Cross-validation statistics can then be analysed with other methods
areModelsDominated(summary(runs))
cast(melt(summary(runs),id.vars=c("calib.period","sim.period")),calib.period~variable+sim.period)
paretoTimeAnalysis(runs)
}
