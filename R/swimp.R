## SWIMP: Simple Wetland Inundation Model using Poweroids
## by Felix Andrews, Barry Croke and Baihua Fu 2009
## (inspired by MFAT Floodplain Hydrology model)



#' Simple Wetland Inundation Model using Poweroids
#'
#' Model flood area / duration / depth in wetlands.
#' %% ~~ If necessary, more details than the description above ~~
#'
#'
#'
#' @param flow.ML inflow or streamflow in ML per timestep.
#' @param thresh a threshold for \code{flow.ML}, such that only flow above this
#' value enters the wetland.
#' @param alpha,beta parameters defining the shape of the wetland. See
#' \code{\link{poweroid}}.
#' @param E.mm,P.mm potential evapo-transpiration and precipitation in mm per
#' timestep.
#' @param Ksat.mm.day Saturated hydraulic conductivity in mm per timestep,
#' relative to a reference pressure of 10cm. If this is 0, the wetland surface
#' water is isolated from the surrounding water table, i.e. there is no
#' infiltration nor discharge.
#' @param e Placeholder
#' @param g stress threshold in terms of Catchment Moisture Deficit (mm), as in
#' the IHACRES CMD model, where \code{g = f * d}. See
#' \code{\link{IHACRES.CMD.model}}.
#' @param Hmax Placeholder
#' @param Amax Placeholder
#' @param porosity effective porosity of the soil.
#' @param M_0 initial value of Catchment Moisture Deficit, mm.
#' @param V_0 initial volume of surface water in wetland, ML.
#' @param drainage drainage rate as a proportion of volume above
#' \code{drainLevel} per timestep.
#' @param drainLevel water level (millimetres from base) above which drainage
#' occurs.
#' @return a \code{\link{zoo}} object (time series object).
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{poweroid}}, \code{\link{convertFlow}}
#' @references ...
#' @keywords math models
#' @examples
#'
#' ## assume Q is inflow in ML/day
#' set.seed(1)
#' Q <- rpois(100, lambda = 0.1) * 1000
#' ## assume depth distribution follows a cone, i.e. beta = 1
#' ## estimate alpha given known area vs volume
#' ## lets say a volume of 1000 ML corresponds to area 20 km^2
#' alpha <- poweroid(V = 1000, A = 20, beta = 1)$alpha
#' flood <- swimp(Q, alpha = alpha, beta = 1, E.mm = 10)
#' head(flood, 20)
#' xyplot(flood)
#' @useDynLib hydromad swimp_core
#' @export
swimp <-
  function(flow.ML,
           thresh = 0,
           alpha,
           beta,
           E.mm = 0,
           P.mm = 0,
           Ksat.mm.day = 0,
           e = 0.2,
           g = 140,
           Hmax = 2000,
           Amax = 10000,
           porosity = 0.2,
           M_0 = Hmax * porosity,
           V_0 = 0,
           drainage = 0,
           drainLevel = 0) {
    stopifnot(NROW(flow.ML) > 1)
    stopifnot(NCOL(flow.ML) == 1)
    force(list(alpha, beta))
    overflow <- pmax(coredata(flow.ML) - thresh, 0)
    E.mm <- rep(coredata(E.mm), length = length(overflow))
    P.mm <- rep(coredata(P.mm), length = length(overflow))
    ## volume of surface water, ML
    V <- double(length(overflow))
    ## area of surface water, km^2
    A <- double(length(overflow))
    ## surface water level, mm
    H <- double(length(overflow))
    ## ground water level as moisture deficit (CMD), mm
    M <- double(length(overflow))
    ## infiltration into wetland
    Iw <- double(length(overflow))
    ## initial conditions
    M[1] <- M_0
    V[1] <- V_0
    ## work out volume below drainage level
    undrainedV <- poweroid(H = drainLevel, alpha = alpha, beta = beta)$V

    ## TODO: start from t=1 (add dummy row to time series?)
    if (isTRUE(hydromad.getOption("pure.code"))) {
      ## slower version in R for cross-checking
      storedV <- poweroid(H = drainLevel, alpha = alpha, beta = beta)$V
      for (t in 2:NROW(overflow)) {
        ## work out water level in ground
        ## relative to base of wetland, H==0
        Hg <- Hmax - (1 / porosity) * M[t - 1]
        ## infiltration rate (drainage) from wetland
        ## depends on pressure; reference = 100mm
        Iw[t] <- Ksat.mm.day * (H[t - 1] - Hg) / 100
        ## if discharging into wetland, adjust mass for porosity
        if (H[t - 1] < Hg) Iw[t] <- Iw[t] * porosity
        ## evapo-transpiration (equation from CMD model)
        Mf <- M[t - 1] - P.mm[t]
        ETg <- e * E.mm[t] * min(1, exp(2 * (1 - Mf / g)))
        ETg <- max(ETg, 0)
        ## mass balance of water level in ground (CMD)
        M[t] <- M[t - 1] + ETg - P.mm[t] - Iw[t] * (A[t - 1] / Amax)
        M[t] <- max(M[t], 0)
        ## mass balance of volume in wetland
        V[t] <- V[t - 1] + overflow[t] + (P.mm[t] - E.mm[t]) * A[t - 1] - Iw[t] * A[t - 1]
        V[t] <- max(V[t], 0)
        ## drainage of volume above drainage level
        V[t] <- V[t] - max(V[t] - undrainedV, 0) * drainage
        ## convert volume to water level and area
        ## NOTE: these formulae are designed to avoid numerical explosions!
        ## e.g. ^(1/beta) when beta = 0.01
        H[t] <- ((2 / beta + 1) * V[t] / pi)^(beta / (2 + beta)) * alpha^(2 / (2 + beta))
        A[t] <- ((2 / beta + 1) * V[t] / (alpha * pi^(-beta / 2)))^(2 / (2 + beta))
      }
    } else {
      ans <- .C(swimp_core,
        as.double(overflow),
        as.integer(length(overflow)),
        as.double(alpha),
        as.double(beta),
        as.double(E.mm),
        as.double(P.mm),
        as.double(Ksat.mm.day),
        as.double(e),
        as.double(g),
        as.double(Hmax),
        as.double(Amax),
        as.double(porosity),
        as.double(drainage),
        as.double(undrainedV),
        V = V,
        A = A,
        H = H,
        M = M,
        Iw = Iw,
        PACKAGE = "hydromad"
      )
      V <- ans$V
      A <- ans$A
      H <- ans$H
      M <- ans$M
      Iw <- ans$Iw
    }
    mean.depth <- ifelse(A > 0, V / A, 0)
    ans <- zoo(
      cbind(
        volume = V, area = A,
        level = H, mean.depth = mean.depth
      ),
      time(flow.ML)
    )
    if (Ksat.mm.day > 0) {
      ans <- cbind(ans, Iw = Iw, CMD = M)
    }
    ans
  }

#' @import utils
utils::globalVariables(c("swimp_core"))
