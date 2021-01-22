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
#' T (average air temperature in ºC) and E (potential
#' evapotranspiration in mm). If E is not supplied then the PET argument must
#' be used.
#' @param PET A list containing named objects "PET" and optionally "Tmean".
#' PET and Tmean must be vectors of length 12 or 365 that represent mean values.
#' Optional if `E` is supplied in `DATA`. If supplied then `cet` must be
#' parameterised.
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
#' @param cet Potential ET correction factor. Optional if a full PET series
#' is provided.
#' @param initialise_sm If true, the soil moisture store is initialised to
#' equal `fc`*`lp` to match HBVlight behaviour. Defaults to false.
#' @param return_state Whether to return the state variables.
#'
#' @param U Effective rainfall series/recharge series.
#' @param perc Maximum percolation from upper to lower groundwater storage.
#' @param uzl Threshold for quick runoff for k0 outflow (mm).
#' @param k0 Recession coefficient (quick runoff).
#' @param k1 Recession coefficient (upper groundwater storage).
#' @param k2 Recession coefficient (lower groundwater storage).
#' @param maxbas Routing, length of triangular weighting function (days).
#' @param epsilon Values smaller than this in the output will be set to zero.
#' @param initial_slz Initial value for the lower store (SLZ). Defaults to 0.
#' @param return_components Whether to return state variables from the routing
#' routine.
#' @details This implementation of this HBV model closely follows the
#' description of HBV light by Seibert and Vis, 2012. Daily average temperature
#' data is required for the snow routine. If daily potential evapotranspiration
#' (PET) data is not provided in `DATA`, then the PET is estimated using the
#' HBV method and the `PET` and `cet` arguments must be specified. The list
#' needs to contain a vector named `"PET"` containing daily average (or matching
#' the timestep of `DATA`) values of length 12 (monthly) or 365 (days of year).
#' `"Tmean"` can also be provided in this list, of length 12 or 365, otherwise
#' average monthly values will be calculated from the average temperature
#' series in `DATA`. See an example below of how to pass average values.
#' @return The timeseries of simulated streamflow (U). If return state is set
#' to true, the state variables of the model are also returned. These include:
#' snow depth (Snow), soil moisture (SM), potential evapotranspiration (PET)
#' and actual evapotranspiration (AET), .
#'
#' For hbv_routing, the routed effective rainfall series (X) is returned. If
#' return components is to true, the state variables of the routing model are
#' returned: upper groundwater storage (SUZ), lower groundwater storage (SLZ),
#' runoff from quick flow (Q0), upper groundwater (Q1) and lower groundwater
#' (Q2).
#'
#' Default parameter ranges are guided by Seibert (1997) and Seibert and Vis
#' (2012) and (see the references section). Parameter ranges for your
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
#'
#' # If using average values of PET and Tmean if daily values of PET are not
#' # available and E is not in DATA. cet must be specified.
#' \dontrun{
#' mod <- hydromad(
#'   DATA = Corin,
#'   sma = "hbv",
#'   routing = "hbvrouting",
#'   PET = list("PET" = PET, "Tmean" = Tmean),
#'   cet = 0.1
#' )
#' }
#'
#' @keywords models
#' @useDynLib hydromad
#' @importFrom Rcpp sourceCpp
#' @useDynLib hydromad _hydromad_hbv_sim
#' @export

hbv.sim <- function(DATA,
                    tt, cfmax, sfcf, cfr, cwh,
                    fc, lp, beta, cet,
                    return_state = FALSE,
                    initialise_sm = FALSE,
                    PET) {
  # DATA: zoo series with P, Q and T and optionally E
  # PET: list with names PET and optionally Tmean with corresponding vectors
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
  stopifnot(c("P", "T") %in% colnames(DATA))

  # If daily PET isn't provided calculate it
  if (!"E" %in% colnames(DATA)) {
    stopifnot("PET" %in% names(PET))
    if (missing(cet)) stop("missing parameter cet")
    if ("Tmean" %in% names(PET)) {
      DATA <- hbv.pet(
        DATA = DATA,
        PET = PET$PET,
        Tmean = PET$Tmean,
        cet = cet
      )
    } else {
      DATA <- hbv.pet(
        DATA = DATA,
        PET = PET$PET,
        cet = cet
      )
    }
  }
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
  stopifnot(is.logical(initialise_sm))

  inAttr <- attributes(DATA[, 1])
  DATA <- as.ts(DATA)

  P <- DATA[, "P"]
  Tavg <- DATA[, "T"]
  E <- DATA[, "E"]

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
    E[bad] <- NA
    AET[bad] <- NA

    attributes(sp) <- inAttr
    attributes(sm) <- inAttr
    attributes(E) <- inAttr
    attributes(AET) <- inAttr

    ans <- cbind(
      U = U,
      Snow = sp,
      SM = sm,
      PET = E,
      AET = AET
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
                           perc, uzl,
                           k0, k1, k2,
                           maxbas,
                           initial_slz = I(0),
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
      wi, n_maxbas, initial_slz
    )
    suz <- ans$suz
    slz <- ans$slz
    Q0 <- ans$Q0
    Q1 <- ans$Q1
    Q2 <- ans$Q2
  } else { # R version
    # Initialise variables and vectors
    suz_ <- 0
    slz_ <- initial_slz

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
  Q0 <- zoo::rollapplyr(
    c(rep(0, n_maxbas - 1), Q0), n_maxbas, function(Q) sum(Q * wi)
  )
  Q1 <- zoo::rollapplyr(
    c(rep(0, n_maxbas - 1), Q1), n_maxbas, function(Q) sum(Q * wi)
  )
  Q2 <- zoo::rollapplyr(
    c(rep(0, n_maxbas - 1), Q2), n_maxbas, function(Q) sum(Q * wi)
  )

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

# Internal function to calculate potential ET using average PET and
# and temperature
#' @useDynLib hydromad _hydromad_hbv_pet
hbv.pet <- function(DATA, PET, Tmean, cet) {
  # DATA: the DATA series that goes to hbv.sim
  # PET: mean PET. Vector of length 12 or 365
  # Tmean: mean temperature. Vector of length 12 or 365. If left out
  # average temp is calculated from DATA$T
  # cet: parameter to adjust PET
  stopifnot(length(PET) %in% c(12, 365))
  stopifnot("T" %in% colnames(DATA))

  dates <- index(DATA)
  Tavg <- DATA$`T`

  if (missing(Tmean)) {
    # If Tmean isn't provided then calculate mean monthly temperature
    mon <- as.integer(strftime(dates, "%m"))
    Tmean <- as.vector(aggregate(Tavg, by = mon, mean, na.rm = TRUE))
  } else {
    stopifnot(length(Tmean) %in% c(12, 365))
  }

  if (length(Tmean) == 12) {
    # idx correspond to month midpoints From Dec -> Jan-Dec -> Jan
    idx <- c(1, 32, 62, 91, 122, 152, 183, 213, 244, 275, 305, 336, 366, 397)
    Tmean <- stats::approxfun(idx, y = c(Tmean[12], Tmean, Tmean[1]))(17:381)
  }

  if (length(PET) == 12) {
    idx <- c(1, 32, 62, 91, 122, 152, 183, 213, 244, 275, 305, 336, 366, 397)
    PET <- stats::approxfun(idx, y = c(PET[12], PET, PET[1]))(17:381)
  }

  COMPILED <- hydromad.getOption("pure.R.code") == FALSE
  if (COMPILED) {
    E <- hbv_pet(
      dates = dates,
      Tavg = Tavg,
      PET = PET,
      Tmean = Tmean,
      cet = cet
    )
  } else {
    E <- rep(0, length(Tavg))
    # Deal with leap years. Feb 29 = Feb 28
    doy <- as.integer(strftime(dates, "%j"))
    year <- as.integer(strftime(dates, "%Y"))
    idx <- ((year %% 4 == 0) & ((year %% 100 != 0) | (year %% 400 == 0))) &
      doy > 59
    doy[idx] <- doy[idx] - 1
    for (t in 1:365) {
      # Mean temperature for doy
      # Eq7 in Seibert and Vis 2012
      idx <- which(doy == t)
      pet_ <- 1 + cet * (Tavg[idx] - Tmean[t])
      pet_[pet_ > 2] <- 2
      pet_[pet_ < 0] <- 0
      E[idx] <- pet_ * PET[t]
    }
  }
  DATA$E <- E
  return(DATA)
}
