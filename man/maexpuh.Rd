\name{maexpuh}
\alias{maexpuh}
\alias{maexpuh.sim}
\title{
  Smoothed exponential stores
}
\description{
  Multiple exponential stores with pre- and post- smoothing on 2nd and 3rd stores
}
\usage{
maexpuh.sim(U, delay = 0, tau_s = 0, tau_q = 0, tau_3 = 0, v_s = 1, 
    v_q = NA, v_3 = 0, series = 0, loss = 0, Xs_0 = 0, Xq_0 = 0, 
    X3_0 = 0, w_s = 7, w_3 = 30, pars = NULL, return_components = FALSE, 
    na.action = na.pass, epsilon = hydromad.getOption("sim.epsilon"))
}
\arguments{
  \item{U}{
}
  \item{delay}{
}
  \item{tau_s}{
}
  \item{tau_q}{
}
  \item{tau_3}{
}
  \item{v_s}{
}
  \item{v_q}{
}
  \item{v_3}{
}
  \item{series}{
}
  \item{loss}{
}
  \item{Xs_0}{
}
  \item{Xq_0}{
}
  \item{X3_0}{
}
  \item{w_s}{
}
  \item{w_3}{
}
  \item{pars}{
}
  \item{return_components}{
}
  \item{na.action}{
}
  \item{epsilon}{
}
}
\details{
}
\value{
}
\author{
  Joseph Guillaume
}
\seealso{
  \code{\link{expuh}}
}
\examples{
  ## Pre and post filtered Xs
  Us[]=filter(U,rep(1/w_s,w_s),sides=1)
  Us[1:(w_s-1)] <- cumsum(U[1:(w_s-1)])/1:(w_s-1)
  Xstemp[] <- filter_loss(beta_s * U, alpha_s, loss = lossVal,
                      init = Xs_0)
  Xs[]=filter(Xstemp,rep(1/w_s,w_s),sides=1)
  Xs[1:(w_s-1)] <- cumsum(Xstemp[1:(w_s-1)])/1:(w_s-1)
}
\keyword{ ts }

