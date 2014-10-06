\name{paretoObjectivesNsga2}
\alias{paretoObjectivesNsga2}
\title{
  Multi-objective optimisation using NSGAII
}
\description{
  Estimate multi-objective Pareto front using NSGAII
}
\usage{
paretoObjectivesNsga2(MODEL, objective = hydromad.getOption("objective"), control = hydromad.getOption("nsga2.control"))
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
  \item{control}{
    arguments for nsga2 function. See \code{\link{nsga2}}
}
}
\value{
  \code{\link{runlist}} of models on Pareto front
}
\author{
Joseph Guillaume
}
\examples{
data(Cotter)
x <- Cotter[1:1000]

## IHACRES CWI model with exponential unit hydrograph
## an unfitted model, with ranges of possible parameter values
modx <- hydromad(x, sma = "cwi", routing = "expuh",
                 tau_s = c(2,100), v_s = c(0,1))
## Multi-objective optimisation
front <- paretoObjectivesNsga2(modx,objective=list(hmadstat("r.sq.log"),hmadstat("r.squared")))
## Pairwise plot of parameters on Pareto front
splom(coef(front))
## Calculate objectives
stats <- t(sapply(front,objFunVal,objective=list(hmadstat("r.sq.log"),hmadstat("r.squared"))))
splom(stats)
}
\keyword{optimization }