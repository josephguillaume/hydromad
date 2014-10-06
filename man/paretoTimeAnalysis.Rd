\name{paretoTimeAnalysis}
\alias{paretoTimeAnalysis}
\title{
Pareto analysis of performance across time periods
}
\description{
  Produces summary results assessing how performance of non-dominated
  models varies across periods
}
\usage{
paretoTimeAnalysis(rl, objectives = "r.squared", qoi = list(q90 = function(Q, X, ...) quantile(X, 0.9, na.rm = TRUE), r.sq.log = hmadstat("r.sq.log")))
}
\arguments{
  \item{rl}{
    a crossvalidation \code{runlist}: a list of fitted model objects
    produced by \code{\link{crossValidate}}
}
  \item{objectives}{
    Vector of column names containing performance measures to be used by
    \code{\link{areModelsDominated}} to
    determine whether models are dominated across time periods.
}
  \item{qoi}{
    Quantities of interest to calculate, as interpreted by
    \code{\link{objFunVal}}. Defaults are given as examples: 90th
    percentile flow (a scalar prediction) and R Squared using log transformed
    data (a performance statistic).
}
}
\value{
  None, used for its side effect of producing text.
}
\references{
}
\author{
Joseph Guillaume
}
\seealso{
  \code{\link{crossValidate}} for an example use, 
  \code{\link{paretoCatchments}} for the first part of the analysis
  using a data.frame of performance statistics
}