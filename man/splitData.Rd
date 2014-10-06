\name{splitData}
\alias{splitData}
\title{
  Split data of a model instance
}
\description{
  Creates a runlist of model instances with separate data periods
}
\usage{
splitData(object, start.dates=NULL, end.dates=NULL,periods=NULL)
}
\arguments{
  \item{object}{
    an object of class \code{hydromad}.
}
  \item{start.dates}{
    dates at which each model instance should start, passed as \code{start} to \code{\link{window.zoo}}
}
\item{end.dates}{
      dates at which each model instance should end, passed as \code{end} to \code{\link{window.zoo}}
    }
    \item{periods}{
      named list of start and end dates. Overrides \code{start.dates} and
      \code{end.dates} if it is not \code{NULL}
    }
}
\value{
  A runlist of hydromad objects with data periods starting at
  \code{start.dates} and ending at \code{end.dates}
}
\author{
  Joseph Guillaume
}
\seealso{
  \code{\link{crossValidate}}
}
\examples{
data(Cotter)
modx <- hydromad(Cotter, sma = "cwi", routing = "expuh",
                 tau_s = c(2,100), v_s = c(0,1))
start.dates=as.Date(c("1966-05-01","1976-05-1"))
end.dates=as.Date(c("1976-04-30","1986-04-30"))
rl=splitData(modx,start.dates,end.dates)
str(rl)
## Fit the first period
fitx <- fitByOptim(rl[[1]])
## Evaluate the resulting model on all periods
rl.first=update(rl,newpars=coef(fitx))
}