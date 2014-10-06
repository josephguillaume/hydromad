\name{paretoObjectivesVaryWeights}
\alias{paretoObjectivesVaryWeights}
\title{
  Multi-objective optimisation by varying weights
}
\description{
  Estimate multi-objective Pareto front using multiple weighted single
  objective optimisations
}
\usage{
paretoObjectivesVaryWeights(MODEL, objective = hydromad.getOption("objective"), weights, fitBy, ...)
}
\arguments{
  \item{MODEL}{
        a model specification created by \code{\link{hydromad}}. It should
    not be fully specified, i.e one or more parameters should be defined
    by \emph{ranges} of values rather than exact values.
}
\item{objective}{
      objective functions to maximise, as a list with elements as \code{function(Q, X, ...)}.
    See \code{\link{objFunVal}}.
}
\item{weights}{
  matrix of weights to use, with columns in the same order as \code{objective}
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
  \code{\link{runlist}} of models on Pareto front, corresponding to the
  \code{weights} specified for each \code{objective}
}
\author{
Joseph Guillaume
}
\seealso{
  \code{\link{paretoObjectivesNsga2}}
}
\examples{
data(Cotter)
x <- Cotter[1:1000]

## IHACRES CWI model with exponential unit hydrograph
## an unfitted model, with ranges of possible parameter values
modx <- hydromad(x, sma = "cwi", routing = "expuh",
                 tau_s = c(2,100), v_s = c(0,1))

## Uncomment to parallelise the fitBy runs
#library(parallel)
#cl <- makeCluster(getOption("cl.cores", 2))
#clusterEvalQ(cl,library(hydromad))
#hydromad.options(parallel="clusterApply")


## Optimisation of multiple weights for r.sq.log and r.squared
weights<-cbind(c(0,0.33,0.5,0.67,1),1-c(0,0.33,0.5,0.67,1))

## Estimate parameters using single fitByOptim run
## from single initial parameter set
front <- paretoObjectivesVaryWeights(modx,objective=list(hmadstat("r.sq.log"),hmadstat("r.squared")),
  weights=weights,fitBy=fitByOptim,samples=1)

summary(front)
## Plot objectives
stats <- t(sapply(front,objFunVal,objective=list(hmadstat("r.sq.log"),hmadstat("r.squared"))))
plot(stats)
}
\keyword{optimization }