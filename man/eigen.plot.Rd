% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eigenplot.R
\name{eigen.plot}
\alias{eigen.plot}
\alias{eigen.plot.single}
\title{Eigenplot}
\usage{
eigen.plot(obj, fixed.axis = TRUE)
}
\arguments{
\item{obj}{\code{rsm} object, as produced by \code{\link{evalRSM}}}

\item{fixed.axis}{Whether to use the same dimensions for each pair of
variables. Allows comparison of relative width/height of ellipses}
}
\value{
Plot of ellipses for each pair of variables, with larger
width/height of ellipse indicating lower level of identifiability, and
rotation of ellipse away from the vertical indicating level of interaction.
}
\description{
Plot pair-wise eigen values of quadratic response surface model
}
\seealso{
\code{\link{evalRSM}} which includes an example of \code{eigenplot}
}
\author{
Dario Mavec, Mun-Ju Shin, Joseph Guillaume
}
\keyword{models}
