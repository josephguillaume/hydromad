## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Joseph Guillaume <joseph.guillaume@anu.edu.au>
##

#' Exponential components transfer function models with layered slowflow stores
#'
#' A unit hydrograph with a quickflow pathway and two layered slowflow pathways
#' modelling recharge to groundwater in order to allow modelling of long-term
#' disconnection of slowflow stores from streamflow.
#'
#'
#' The \code{expuh3s} model consists of a single quickflow pathway modelled as
#' an exponential store, and a slowflow pathway comprised of two layered
#' stores.
#'
#' Each slowflow store is modelled as a \code{\link{leakyExpStore}}, which has
#' a loss term, has no flow when the store drops below a given level, and can
#' therefore model longer-term disconnection of a store from streamflow.
#'
#' Adapted from Herron and Croke (2009):
#'
#' The upper store, G1, receives rainfall inputs and discharges to the stream,
#' Qs and recharges the lower store. G1 has a lower limit of 0, where flow
#' ceases representing the fully 'drained' condition. Conceptually, the upper
#' store can be viewed as a perched water table, which develops in response to
#' rain and tends to be relatively short-lived, perhaps seasonal. Thus the time
#' constant, \code{tau_s}, for discharge from the 'soil' store will be
#' somewhere between that for quickflow, \code{tau_q} and the groundwater
#' discharge constant, \code{tau_g}.
#'
#' G2 is recharged from G1 when \code{G1>G_1} and discharges to the stream
#' \code{Q_g} when \code{G2>0}. The sum of \code{Q_s} and \code{Q_g} represents
#' the total slowflow pathway. We assume that all extraction and natural
#' groundwater losses (\code{loss}) are from G2. The approach avoids the need
#' to specify a maximum capacity for either storage, but the introduction of a
#' recharge term, \code{R} between the stores adds a new parameter.
#'
#' Recharge is represented by a constant rate \code{R} which ceases when
#' \code{G1<G_1}, diminishing linearly to that point when
#' \code{thres<G1<thres+loss}. Setting \code{G_1=0} (the default) ceases
#' recharge when flow ceases.
#'
#' @name expuh3s.sim
#' @param U input time series (units below assume ML/day)
#' @param delay lag (dead time) between input and response, in time steps.
#' @param v_s Fraction of effective rainfall that goes to groundwater
#' @param tau_q Recession coefficient for quickflow (days)
#' @param tau_s Recession coefficient for soil store (G_1) discharge (days)
#' @param tau_g Recession coefficient for groundwater store (G_2) discharge
#' (days)
#' @param R Maximum recharge from G_1 to G_2 (ML/day)
#' @param G_1 storage threshold to stop recharge (ML) (less than zero)
#' @param loss Groundwater loss (ML/day)
#' @param G_2 storage threshold to stop groundwater loss (ML) (less than zero)
#' @param Xs_0,Xq_0,Xg_0 initial values of the exponential components.
#' @param pars the parameters as a named vector. If this is given, it will
#' over-ride the named parmameter arguments.
#' @param return_components whether to return all component time series.
#' @param na.action function to remove missing values, e.g.
#' \code{\link[=na.omit.ts]{na.omit}}.
#' @param epsilon values smaller than this in the output will be set to zero.
#' @return the model output as a \code{\link{ts}} object, with the same
#' dimensions and time window as the input \code{U}.  If
#' \code{return_components = TRUE}, it will have multiple columns named
#' \code{Xs}, \code{Xq} and \code{Xg}.
#' @author Joseph Guillaume \email{joseph.guillaume@@anu.edu.au}
#' @seealso \code{\link{expuh}},\code{\link{leakyExpStore}}
#' @references Herron, N.F. and B.F.W. Croke (2009). IHACRES-3S - A 3-store
#' formulation for modelling groundwater-surface water interactions. In
#' Anderssen, R.S., R.D. Braddock and L.T.H. Newham (eds) \emph{18th World
#' IMACS Congress and MODSIM09 International Congress on Modelling and
#' Simulation.} Modelling and Simulation Society of Australia and New Zealand
#' and International Association for Mathematics and Computers in Simulation,
#' July 2009, pp. 3081-3087. ISBN: 978-0-9758400-7-8.
#' \url{http://www.mssanz.org.au/modsim09/I1/herron.pdf}
#' @keywords ts
#'
#' @export
expuh3s.sim <-
  function(U, delay = 0, v_s,
           tau_q, tau_s, tau_g,
           R, G_1, loss, G_2,
           Xs_0 = 0, Xq_0 = 0, Xg_0 = 0,
           pars = NULL,
           return_components = FALSE,
           na.action = na.pass,
           epsilon = hydromad.getOption("sim.epsilon")) {
    if (!is.null(pars)) {
      ccall <- match.call()
      ccall$pars <- NULL
      ccall <- as.call(modifyList(
        as.list(ccall),
        as.list(pars)
      ))
      return(eval.parent(ccall))
    }
    delay <- round(delay)
    ## note U is allowed to be multi-variate, i.e. multiple columns
    U <- na.action(U)
    ## apply 'U' delay in reverse to 'X' (i.e. lag by 'delay' steps)
    ## so can take delay as 0 for simulation purposes
    if (delay != 0) {
      U <- lag(U, -delay)
    }

    stopifnot(all(c(tau_s, tau_q, tau_g) >= 0))

    ## convert from 'tau' and 'v' to 'alpha' and 'beta'
    alpha_q <- exp(-1 / tau_q)
    beta_q <- (1 - v_s) * (1 - alpha_q)

    ### apply exponential decay filter to quick flow
    Xq <- U * NA
    Xq[] <- filter_ok(beta_q * U, alpha_q, method = "recursive", init = Xq_0)

    ## Upper store
    upper <- leakyExpStore.sim(v_s * U, tau_s, loss = R, thres = G_1, init = Xs_0, return_components = TRUE)
    Xs <- upper[, "Q"]

    ## Lower store
    lower <- leakyExpStore.sim(upper[, "L"], tau_g,
      loss = loss, thres = G_2, init = Xg_0,
      return_components = FALSE
    )
    Xg <- lower[, "Q"]

    ## align results to original input
    Xs <- shiftWindow(Xs, delay)
    Xq <- shiftWindow(Xq, delay)
    Xg <- shiftWindow(Xg, delay)

    ## zap simulated values smaller than epsilon
    Xs[Xs < epsilon] <- 0
    Xq[Xq < epsilon] <- 0
    Xg[Xg < epsilon] <- 0

    ## TODO: return store as well as flows
    if (return_components) {
      return(cbind(Xs = Xs, Xq = Xq, Xg = Xg))
    } else {
      return(Xs + Xq + Xg)
    }
  }
