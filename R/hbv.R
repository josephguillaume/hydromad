# HBV model code for hydromad
# R and Rcpp code by Alexander Buzacott (abuz5257@uni.sydney.edu.au)
#
# Implementation based on HBV light as described in:
# Seibert, J. and Vis, M. (2012). Teaching hydrological modeling with a user-
# friendly catchment runoff-model software package. Hydrology and Earth System
# Sciences, 16, 3315–3325, 2012.
#
# HBV references:
# Bergström, S. and Forsman, A.: Development of a Conceptual Deterministic
# Rainfall-Runoff Model, Nordic Hydrology, 4(3), 147–170, 1973.
#
# Bergström, S.: The HBV Model: Its Structure and Applications,Swedish
# Meteorological and Hydrological Institute (SMHI), Hydrology, Norrköping, 35
# pp., 1992.
#
# Bergström, S.: The HBV model (Chapter 13), in: Computer Models of Watershed
# Hydrology, edited by: Singh, V. P., Water Resources Publications, Highlands
# Ranch, Colorado, USA, 443–476, 1995.

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
#' @param cwh Liquid water holding capacity of the snowpack.
#' @param cfr Refreezing coefficient for water in the snowpack.
#' @param fc Maximum amount of soil moisture storage (mm).
#' @param lp Threshold for reduction of evaporation. Limit for potential
#' evapotranspiration.
#' @param beta Shape coefficient in soil routine.
#' @param k0 Recession coefficient (quick runoff).
#' @param k1 Recession coefficient (upper groundwater storage).
#' @param k2 Recession coefficient (lower groundwater storage).
#' @param uzl Threshold for quick runoff for k0 outflow (mm).
#' @param perc Maximum percolation from upper to lower groundwater storage.
#' @param maxbas Routing, length of triangular weighting function.
#' @param return_state Whether to return the state variables.
#' @details This implementation of this HBV model closely follows the
#' description of HBV light by Seibert and Vis, 2009. Daily average temperature
#' data is required for the snow routine. Daily potential evapotranspiration
#' (PET) data is required as the routine to calculate daily PET using long term
#' mean PET and daily temperature is not included.
#' @details The timeseries of simulated streamflow (U). If return state is set
#' to true, the state variables of the model are also returned. These include:
#' snowpack (sp), water content of snowpack (wc), soil moisture (sm), actual
#' evapotranspiration (ETa), upper groundwater storage (uz), lower groundwater
#' storage (lz). Note, there are no other components to return with
#' `hbvrouting`.
#' @references
#' Seibert, J. and Vis, M. (2012). Teaching hydrological modeling with a
#' user-friendly catchment-runoff-model software package. Hydrology and Earth
#' System Sciences, 16, 3315–3325, 2012.
#'
#' Bergström, S.: The HBV Model: Its Structure and Applications, Swedish
#' Meteorological and Hydrological Institute (SMHI), Hydrology, Norrköping,
#' 35 pp., 1992.
#'
#' Bergström, S.: The HBV model (Chapter 13), in: Computer Models of
#' Watershed Hydrology, edited by: Singh, V. P., Water Resources Publications,
#' Highlands Ranch, Colorado, USA, 443–476, 1995.
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
#' @useDynLib hydromad
#' @importFrom Rcpp sourceCpp
#' @useDynLib hydromad _hydromad_hbv_sim
#' @export

hbv.sim <- function(DATA,
                    tt, cfmax, sfcf, cwh, cfr,
                    fc, lp, beta,
                    k0, k1, k2, uzl, perc,
                    return_state = FALSE) {
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

  # Groundwater routine
  # k0: recession coefficient
  # k1: recession coefficient
  # k2: recession coefficient
  # uzl: upper zone layer threshold
  # perc: percolation from upper to lower response box

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
  stopifnot(k0 >= 0)
  stopifnot(k1 >= 0)
  stopifnot(k2 >= 0)
  stopifnot(uzl >= 0)
  stopifnot(perc >= 0)

  inAttr <- attributes(DATA[, 1])
  DATA <- as.ts(DATA)

  P <- DATA[, "P"]
  # Not needed Q <- DATA[, "Q"]
  E <- DATA[, "E"]
  Tavg <- DATA[, "T"]

  # Skip missing values
  bad <- is.na(P) | is.na(E) | is.na(Tavg)
  P[bad] <- 0
  E[bad] <- 0
  Tavg[bad] <- 0

  # Check for C++
  COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
  if (COMPILED) {
    # Run C++
    ans <- hbv_sim(
      P, E, Tavg,
      tt, cfmax, sfcf, cwh, cfr,
      fc, lp, beta,
      k0, k1, k2, uzl, perc
    )
    U <- ans$U
    if (return_state == TRUE) {
      ETa <- ans$ETa
      sm <- ans$sm
      sp <- ans$sp
      wc <- ans$wc
      uz <- ans$uz
      lz <- ans$lz
    }
  } else { # Run R Model
    # Set up vectors
    sm <- rep(0, nrow(DATA)) # Soil water storage
    sp <- rep(0, nrow(DATA)) # Snow store
    wc <- rep(0, nrow(DATA)) # Depth of liquid in snow store
    uz <- rep(0, nrow(DATA)) # Shallow soil storage
    lz <- rep(0, nrow(DATA)) # Deep soil storage
    Qsim <- rep(0, nrow(DATA)) # Flow from reservoirs
    ETa <- rep(0, nrow(DATA)) # Actual ET

    # Run model, starting at day 2
    for (t in 2:nrow(DATA)) {
      # ------------------------------------------------------------------------
      # Snow routine
      # ------------------------------------------------------------------------
      infil <- 0

      # Determine if P is snow or rain
      if (Tavg[t] < tt) {
        # Refreezing of liquid in snow pack
        refr <- min(c(cfr * cfmax * (Tavg[t] - tt), wc[t - 1]))
        wc[t] <- wc[t - 1] - refr
        # Use snowfall correction factor
        snow <- P[t] * sfcf
        # Add snow to snow pack + refreezing water
        sp[t] <- sp[t - 1] + snow + refr
        infil <- P[t] - snow
      } else {
        # Calculate and remove snow melt
        melt <- min(cfmax * (Tavg[t] - tt), sp[t - 1])
        sp[t] <- sp[t - 1] - melt
        # Add water to liquid in snow pack
        wc[t] <- wc[t - 1] + melt
        # Calculate maximum liquid water holding capacity of snow pack
        maxwc <- max(sp[t] * cwh, 0)
        if (wc[t] > maxwc) {
          # Add liquid excess water to effective P
          infil <- P[t] + wc[t] - maxwc
          wc[t] <- maxwc
        } else {
          # Retain liquid water in snow pack and just add P
          infil <- P[t]
        }
      }
      # ------------------------------------------------------------------------
      # Soil routine
      # ------------------------------------------------------------------------
      recharge <- 0
      AET <- 0
      runoff <- 0

      # Calculate current soil wetness
      smw <- max(min((sm[t - 1] / fc)^beta, 1), 0)

      # Calculate recharge and take away from infil
      recharge <- infil * smw
      infil <- infil - recharge

      # Add remaining infiltration to soil
      sm[t] <- sm[t - 1] + infil

      # Check if field capacity is exceeded
      if (sm[t] > fc) {
        runoff <- sm[t] - fc
        sm[t] <- fc
      }

      # Calculate actual ET
      AET <- E[t] * min(sm[t] / (fc * lp), 1)
      if (AET < 0) AET <- 0
      # Remove AET from soil if there is water
      if (sm[t] > AET) {
        ETa[t] <- AET
        sm[t] <- sm[t] - AET
      } else {
        ETa[t] <- sm[t]
        sm[t] <- 0
      }

      # ------------------------------------------------------------------------
      # Discharge
      # ------------------------------------------------------------------------
      Q0 <- 0
      Q1 <- 0
      Q2 <- 0

      # Add runoff and recharge to upper zone of storage
      uz[t] <- uz[t - 1] + runoff + recharge

      # Percolation of of water from upper to lower zone
      actPERC <- min(uz[t], perc)
      uz[t] <- uz[t] - actPERC
      lz[t] <- lz[t - 1] + actPERC

      # Calculate runoff from storage
      Q0 <- k0 * max(uz[t] - uzl, 0)
      uz[t] <- uz[t] - Q0

      Q1 <- k1 * uz[t]
      uz[t] <- uz[t] - Q1

      Q2 <- k2 * lz[t]
      lz[t] <- lz[t] - Q2

      Qsim[t] <- Q0 + Q1 + Q2
    } # R loop done
    U <- Qsim
  } # Close ifelse
  # Model run finished

  # Put back missing values
  U[bad] <- NA

  # Attributes
  attributes(U) <- inAttr
  ans <- U

  if (return_state == TRUE) {
    # Return state variables
    sp[bad] <- NA
    wc[bad] <- NA
    sm[bad] <- NA
    ETa[bad] <- NA
    uz[bad] <- NA
    lz[bad] <- NA

    attributes(sp) <- inAttr
    attributes(wc) <- inAttr
    attributes(sm) <- inAttr
    attributes(ETa) <- inAttr
    attributes(uz) <- inAttr
    attributes(lz) <- inAttr

    ans <- cbind(
      U = U,
      ETa = ETa,
      SD = sp,
      LD = wc,
      SM = sm,
      UZ = uz,
      LZ = lz
    )
  }

  return(ans)
}

# Implementation of the triangular weighting function in HBV
# Eq 6 in https://doi.org/10.5194/hess-16-3315-2012
# with a minor adjustment to handle non-integer maxbas

#' @rdname hbv
#' @useDynLib hydromad _hydromad_hbvrouting_sim
#' @export
hbvrouting.sim <- function(U,
                           maxbas,
                           epsilon = hydromad.getOption("sim.epsilon"),
                           return_components = FALSE) {
  # Routing
  # maxbas: routing, length of triangular weighting function

  inAttr <- attributes(U)
  U <- as.ts(U)
  bad <- is.na(U)
  U[bad] <- 0

  ci <- function(u) {
    (2 / maxbas) - abs(u - (maxbas / 2)) * (4 / (maxbas^2))
  }

  n_maxbas <- ceiling(maxbas)
  wi <- rep(0, n_maxbas)

  if (maxbas > 1) {
    for (i in 1:maxbas) {
      wi[i] <- integrate(ci, i - 1, min(i, maxbas))$value[1]
    }
  } else {
    wi <- 1
  }

  COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
  if (COMPILED) {
    X <- hbvrouting_sim(U, wi, n_maxbas)
  } else { # R version
    X <- rep(0, length(U))
    for (t in seq_len(length(X))) { # While t < maxbas, X[t] = U[t]
      if (t < n_maxbas) {
        X[t] <- U[t]
      } else {
        # Else use triangular weighting function
        X[t] <- sum(wi * (U[(t - n_maxbas + 1):t]))
      }
    }
  }

  # Values smaller than epsilon go to 0
  X[abs(X) < epsilon] <- 0
  X[bad] <- NA
  attributes(X) <- inAttr

  return(X)
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
    cwh = c(0, 0.2),
    cfr = c(0, 0.1),
    fc = c(50, 500),
    lp = c(0.3, 1),
    beta = c(1, 6),
    k0 = c(0.05, 0.5),
    k1 = c(0.01, 0.3),
    k2 = c(0.001, 0.1),
    uzl = c(0, 100),
    perc = c(0, 3)
  )
}

#' @rdname hbv
#' @export
hbvrouting.ranges <- function() {
  list(
    maxbas = c(1, 7)
  )
}
