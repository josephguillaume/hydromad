\name{hydromad.options}
\alias{hydromad.options}
\alias{hydromad.getOption}
\title{ User default settings for hydromad }
\description{
  A basic user settings facility, like \code{\link{options}} and
  \code{\link{lattice.options}}.
}
\usage{
hydromad.options(...)
hydromad.getOption(name)
}
\arguments{
  \item{name}{ character giving the name of a setting. }
  \item{...}{ new options can be defined, or existing ones modified, using
          one or more arguments of the form 'name = value' or by
          passing a list of such tagged values.  Existing values can be
          retrieved by supplying the names (as character strings) of
          the components as unnamed arguments.
    }
}
\details{
  These functions are direct copies of the lattice equivalents:
  see \code{\link{lattice.options}}.

  The available options can be seen with \code{str(hydromad.options())}.
  Many of these simply provide defaults for corresponding arguments
  to the \code{\link{hydromad}} function.
}
\seealso{ \code{\link{hydromad}} }
\examples{
oopt <- hydromad.options()
str(oopt)

## reset
hydromad.options(oopt)
}
\keyword{ programming }