# HBV model code for hydromad
# R and Rcpp code by Alexander Buzacott (abuz5257@uni.sydney.edu.au)
#
# Implementation based on HBV light as described in:
# Seibert, J. and Vis, M. (2012). Teaching hydrological modeling with a user-
# friendly catchment runoff-model software package. Hydrology and Earth System
# Sciences, 16, 3315–3325, 2012.
#
#' @name hbv
#' @md
#' @aliases hbv.sim hbvrouting hbvrouting.sim
#' @title HBV rainfall-runoff model
#' @description An implementation of the HBV rainfall-runoff model.
#' @param DATA A time-series like object with columns P (precipitation in mm),
#' E (potential evapotranspiration in mm), T (average air temperature in ºC).
#' @param tt Threshold temperature for snow and snow melt in degrees Celsius.
#' @param cfmax Degree-day factor for snow melt (mm/(ºC.day)).
#' @param sfcf Snowfall correction factor. Amount of precipitation below
#' threshold temperature that should be rainfall instead of snow.
#' @param cfr Refreezing coefficient for water in the snowpack.
#' @param cwh Liquid water holding capacity of the snowpack.
#' @param fc Maximum amount of soil moisture storage (mm).
#' @param lp Threshold for reduction of evaporation. Limit for potential
#' evapotranspiration.
#' @param beta Shape coefficient in soil routine.
#' @param return_state Whether to return the state variables.
#'
#' @param U Effective rainfall series
#' @param perc Maximum percolation from upper to lower groundwater storage.
#' @param uzl Threshold for quick runoff for k0 outflow (mm).
#' @param k0 Recession coefficient (quick runoff).
#' @param k1 Recession coefficient (upper groundwater storage).
#' @param k2 Recession coefficient (lower groundwater storage).
#' @param maxbas Routing, length of triangular weighting function.
#' @param epsilon Values smaller than this in the output will be set to zero.
#' @param return_components Whether to return state variables of the routing
#' routine.
#' @details This implementation of this HBV model closely follows the
#' description of HBV light by Seibert and Vis, 2009. Daily average temperature
#' data is required for the snow routine. Daily potential evapotranspiration
#' (PET) data is required as the routine to calculate daily PET using long term
#' mean PET and daily temperature is not included.
#' @details The timeseries of simulated streamflow (U). If return state is set
#' to true, the state variables of the model are also returned. These include:
#' snow depth (Snow), actual evapotranspiration (AET) , soil moisture (SM).
#'
#' For hbv_routing, if return components is to true, the state variables of
#' the routing model are returned: upper groundwater storage (SUZ), lower
#' groundwater storage (SLZ).
#'
#' Default parameter ranges are guided by Seibert (1997) and Seibert and Vis
#' (2012) and (see the references section), however parameter ranges for your
#' catchment may require either a more restricted or wider range.
#'
#' @references
#'
#' Bergström, S. and Forsman, A.: Development of a Conceptual Deterministic
#' Rainfall-Runoff Model, Nordic Hydrology, 4(3), 147–170, 1973.
#'
#' Bergström, S.: The HBV Model: Its Structure and Applications,Swedish
#' Meteorological and Hydrological Institute (SMHI), Hydrology, Norrköping, 35
#' pp., 1992.
#'
#' Seibert, J. (1997). Estimation of Parameter Uncertainty in the HBV Model.
#' Hydrology Research, 28(4–5), 247–262.
#'
#' Seibert, J. and Vis, M. (2012). Teaching hydrological modeling with a
#' user-friendly catchment-runoff-model software package. Hydrology and Earth
#' System Sciences, 16, 3315–3325, 2012.
#'
#' @author Alexander Buzacott (abuz5257@uni.sydney.edu.au)
#' @seealso `hydromad(sma='hbv', routing='hbvrouting')` to work with
#' models as objects (recommended).
#' @examples
#' # Using example dataset Corin with daily P, Q, potential ET and average T
#' data(Corin)
#'
#' # See default par ranges with hbv.ranges() or hydromad.getOption('hbv')
#' hydromad.getOption("hbv")
#' hydromad.getOption("hbvrouting")
#'
#' # Create model
#' mod <- hydromad(
#'   DATA = Corin,
#'   sma = "hbv",
#'   routing = "hbvrouting"
#' )
#'
#' # Fit using the optim routine with the KGE objective function
#' fit <- fitByOptim(mod, objective = hmadstat("KGE"))
#'
#' # Summary statistics and plot of the fit
#' summary(fit)
#' objFunVal(fit)
#' xyplot(fit)
#' @keywords models
#' @export

hbv.sim <- function(DATA,
                    tt, cfmax, sfcf, cfr, cwh,
                    fc, lp, beta,
                    return_state = FALSE,
                    initialise_sm = FALSE) {
  # DATA: zoo series with P, Q, E and T
  # Snow routine
  # tt: temperature limit for rain/snow (ºC)
  # sfcf: snowfall correction factor
  # cfmax: degree day factor, rate of snow melt (mm/(ºC-d))
  # cfr:  Refreezing factor
  # cwh: Water holding capacity of snow pack

  # Soil routine
  # fc: maximum soil moisture content (mm)
  # lp: limit for potential evapotranspiration
  # beta: parameter in soil routine

  # Check DATA has been entered
  stopifnot(c("P", "E", "T") %in% colnames(DATA))

  # Check valid parameter values have been entered
  stopifnot(is.double(tt))
  stopifnot(is.double(sfcf))
  stopifnot(cfmax >= 0)
  stopifnot(cfr >= 0)
  stopifnot(cwh >= 0)
  stopifnot(fc >= 0)
  stopifnot(lp >= 0)
  stopifnot(beta >= 0)
  stopifnot(cwh >= 0)

  inAttr <- attributes(DATA[, 1])
  DATA <- as.ts(DATA)

  P <- DATA[, "P"]
  # Q <- DATA[, "Q"]
  E <- DATA[, "E"]
  Tavg <- DATA[, "T"]

  # Skip missing values
  bad <- is.na(P) | is.na(E) | is.na(Tavg)
  P[bad] <- 0
  E[bad] <- 0
  Tavg[bad] <- 0

  # Check for C++
  COMPILED <- hydromad.getOption("pure.R.code") == FALSE
  if (COMPILED) {
    # Run C++
    ans <- hbv_sim(
      P, E, Tavg,
      tt, cfmax, sfcf, cfr, cwh,
      fc, lp, beta, initialise_sm
    )
    U <- ans$U
    if (return_state == TRUE) {
      AET <- ans$AET
      sm <- ans$sm
      sp <- ans$sp
    }
  } else { # Run R Model
    # Set up vectors
    sm <- rep(0, nrow(DATA)) # Soil water storage
    sp <- rep(0, nrow(DATA)) # Snow store
    infil <- rep(0, nrow(DATA)) # infiltration to soil
    AET <- rep(0, nrow(DATA)) # Actual ET
    recharge <- rep(0, nrow(DATA)) # Effective precipitation -> water to routing

    # Set up variables
    refr <- 0
    wc_ <- 0
    sp_ <- 0
    sm_ <- ifelse(initialise_sm, fc * lp, 0)

    # Run model
    for (t in seq(1, nrow(DATA))) {
      # ------------------------------------------------------------------------
      # Snow routine
      # ------------------------------------------------------------------------
      infil_ <- 0
      sp_tm1 <- sp_
      # Determine if snow or rain falls
      if (P[t] > 0) {
        if (Tavg[t] > tt) {
          # Precipitation gets added to wc store
          wc_ <- wc_ + P[t]
        } else {
          # Snow and apply snowfall correction factor
          sp_ <- sp_ + P[t] * sfcf
        }
      }
      if (Tavg[t] > tt) {
        # Melt snow
        melt <- cfmax * (Tavg[t] - tt)
        # If melt is greater than snow depth
        if (melt > sp_) {
          # All water is added to infiltration
          infil_ <- sp_ + wc_
          wc_ <- 0
          sp_ <- 0
        } else {
          # Remove melt from snow pack
          sp_ <- sp_ - melt
          wc_ <- wc_ + melt
          # Calculate maximum liquid water holding capacity of snow pack
          maxwc <- sp_ * cwh
          if (wc_ > maxwc) {
            infil_ <- wc_ - maxwc
            wc_ <- maxwc
          }
        }
      } else {
        # Refreeze water in liquid snow store
        refr <- min(cfr * cfmax * (tt - Tavg[t]), wc_)
        sp_ <- sp_ + refr
        wc_ <- wc_ - refr
      }
      sp[t] <- sp_ + wc_
      # infil[t] <- infil_

      # ------------------------------------------------------------------------
      # Soil routine
      # ------------------------------------------------------------------------
      # Divide portion of infiltration that goes to soil/gw
      sm_tm1 <- sm_
      if (infil_ > 0) {
        if (infil_ < 1) {
          infil_s <- infil_
        } else {
          infil_r <- round(infil_)
          infil_s <- infil_ - infil_r
          i <- 1
          while (i <= infil_r) {
            rm <- (sm_ / fc)^beta
            if (rm > 1) rm <- 1
            sm_ <- sm_ + 1 - rm
            recharge[t] <- recharge[t] + rm
            i <- i + 1
          }
        }
        rm <- (sm_ / fc)^beta
        if (rm > 1) rm <- 1
        sm_ <- sm_ + (1 - rm) * infil_s
        recharge[t] <- recharge[t] + rm * infil_s
      }
      # Only AET if there is no snow cover
      if (sp_tm1 == 0) {
        sm_et <- (sm_ + sm_tm1) / 2
        # Calculate actual ET
        AET[t] <- E[t] * min(sm_et / (fc * lp), 1)
        if (AET[t] < 0) AET[t] <- 0
        # Remove AET from soil if there is water
        if (sm_ > AET[t]) {
          sm_ <- sm_ - AET[t]
        } else {
          AET[t] <- sm_
          sm_ <- 0
        }
      }
      sm[t] <- sm_
    } # R loop done
    U <- recharge
  }
  # Put back missing values
  U[bad] <- NA

  # Attributes
  attributes(U) <- inAttr
  ans <- U

  if (return_state == TRUE) {
    # Return state variables
    sp[bad] <- NA
    sm[bad] <- NA
    AET[bad] <- NA

    attributes(sp) <- inAttr
    attributes(sm) <- inAttr
    attributes(AET) <- inAttr

    ans <- cbind(
      U = U,
      Snow = sp,
      SM = sm,
      AET = AET
    )
  }
  return(ans)
}

# Implementation of the triangular weighting function in HBV
# Eq 6 in https://doi.org/10.5194/hess-16-3315-2012
# with a minor adjustment to handle non-integer maxbas

#' @rdname hbv
#' @export
hbvrouting.sim <- function(U,
                           k0, k1, k2, uzl, perc,
                           maxbas,
                           initial_lz = 0,
                           epsilon = hydromad.getOption("sim.epsilon"),
                           return_components = FALSE) {
  # U: effective rainfall series
  # Groundwater routine
  # k0: recession coefficient
  # k1: recession coefficient
  # k2: recession coefficient
  # uzl: upper zone layer threshold
  # perc: percolation from upper to lower response box
  # Routing
  # maxbas: routing, length of triangular weighting function

  stopifnot(k0 >= 0)
  stopifnot(k1 >= 0)
  stopifnot(k2 >= 0)
  stopifnot(uzl >= 0)
  stopifnot(perc >= 0)
  stopifnot(maxbas >= 1)

  inAttr <- attributes(U)
  U <- as.ts(U)
  bad <- is.na(U)
  U[bad] <- 0

  # Calculate maxbas weights
  ci <- function(u) {
    (2 / maxbas) - abs(u - (maxbas / 2)) * (4 / (maxbas^2))
  }

  n_maxbas <- ceiling(maxbas)
  wi <- rep(0, n_maxbas)

  if (maxbas > 1) {
    for (i in 1:n_maxbas) {
      wi[i] <- stats::integrate(ci, i - 1, min(i, maxbas))$value[1]
    }
  } else {
    wi <- 1
  }
  wi <- rev(wi)

  COMPILED <- hydromad.getOption("pure.R.code") == FALSE
  if (COMPILED) {
    ans <- hbvrouting_sim(
      U,
      perc, uzl, k0, k1, k2,
      wi, n_maxbas, initial_lz
    )
    suz <- ans$suz
    slz <- ans$slz
    Q0 <- ans$Q0
    Q1 <- ans$Q1
    Q2 <- ans$Q2
  } else { # R version
    # Initialise variables and vectors
    suz_ <- 0
    slz_ <- initial_lz

    suz <- rep(0, length(U)) # Shallow gw storage
    slz <- rep(0, length(U)) # Deep gw storage

    Q0 <- rep(0, length(U))
    Q1 <- rep(0, length(U))
    Q2 <- rep(0, length(U))

    for (t in seq(1, length(U))) {
      # -----------------------------------------------------------------------
      # Discharge
      # -----------------------------------------------------------------------
      # Add runoff and recharge to upper zone of storage
      suz_ <- suz_ + U[t]
      # Percolation of of water from upper to lower zone
      act_perc <- min(suz_, perc)
      suz_ <- suz_ - act_perc
      slz_ <- slz_ + act_perc

      # Calculate runoff from storage
      Q0[t] <- k0 * max(suz_ - uzl, 0)

      Q1[t] <- k1 * suz_
      suz_ <- suz_ - Q1[t] - Q0[t]

      Q2[t] <- k2 * slz_
      slz_ <- slz_ - Q2[t]

      suz[t] <- suz_
      slz[t] <- slz_
    }
  }
  # Triangular weighting function
  Q0 <- zoo::rollapplyr(Q0, n_maxbas, function(Q) sum(Q * wi), partial = TRUE)
  Q1 <- zoo::rollapplyr(Q1, n_maxbas, function(Q) sum(Q * wi), partial = TRUE)
  Q2 <- zoo::rollapplyr(Q2, n_maxbas, function(Q) sum(Q * wi), partial = TRUE)

  X <- Q0 + Q1 + Q2

  # Values smaller than epsilon go to 0
  X[abs(X) < epsilon] <- 0
  X[bad] <- NA
  attributes(X) <- inAttr

  if (return_components == TRUE) {
    # Return state variables
    suz[bad] <- NA
    slz[bad] <- NA

    Q0[bad] <- NA
    Q1[bad] <- NA
    Q2[bad] <- NA

    attributes(suz) <- inAttr
    attributes(slz) <- inAttr
    attributes(Q0) <- inAttr
    attributes(Q1) <- inAttr
    attributes(Q2) <- inAttr

    ans <- cbind(
      X = X,
      SUZ = suz,
      SLZ = slz,
      Q0 = Q0,
      Q1 = Q1,
      Q2 = Q2
    )
  } else {
    ans <- X
  }
  return(ans)
}

# Suggested parameter ranges
# Guided by https://doi.org/10.2166/nh.1998.15
# and https://doi.org/10.5194/hess-16-3315-2012
#' @rdname hbv
#' @export
hbv.ranges <- function() {
  list(
    tt = c(-2.5, 2.5),
    cfmax = c(1, 10),
    sfcf = c(0.4, 1),
    cfr = c(0, 0.1),
    cwh = c(0, 0.2),
    fc = c(50, 500),
    lp = c(0.3, 1),
    beta = c(1, 6)
  )
}

#' @rdname hbv
#' @export
hbvrouting.ranges <- function() {
  list(
    perc = c(0, 3),
    uzl = c(0, 100),
    k0 = c(0.05, 0.5),
    k1 = c(0.01, 0.3),
    k2 = c(0.001, 0.1),
    maxbas = c(1, 7)
  )
}
