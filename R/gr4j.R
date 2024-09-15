## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

#' GR4J rainfall runoff model
#'
#' GR4J model (modele du Genie Rural a 4 parametres Journalier).
#'
#' The default parameter ranges were taken from the "80% confidence intervals"
#' given in Perrin et. al. (2003).
#'
#' @name gr4j
#' @aliases gr4j.sim gr4jrouting gr4jrouting.sim
#' @param DATA time-series-like object with columns P (precipitation, mm) and E
#' (potential evapo-transpiration, mm).
#' @param U effective rainfall series.
#' @param x1 maximum capacity of the production store (mm).
#' @param x2 groundwater exchange coefficient (mm).
#' @param x3 one day ahead maximum capacity of the routing store (mm).
#' @param x4 time base of unit hydrograph UH1 (time steps).
#' @param etmult Multiplier for the \code{E} input data.
#' @param S_0 Initial soil moisture level as fraction of \code{x1}.
#' @param R_0 Initial groundwater reservoir level as fraction of \code{x3}.
#' @param split Fraction to go into quick flow routing, usually fixed at 0.9.
#' @param return_state to return the series U, S (storage) and ET
#' (evapotranspiration).
#' @param return_components to return the series Xr, Xd and R (reservoir
#' level).
#' @param epsilon values smaller than this in the output will be set to zero.
#' @param transformed transform parameters before use to improve
#' identifiability. They can be untransformed using
#' \code{\link{gr4j.transformpar}}
#' @return the simulated effective rainfall, a time series of the same length
#' as the input series.
#' @author Felix Andrews \email{felix@@nfrac.org} and Joseph Guillaume
#' \email{josephguillaume@@gmail.com}
#' @seealso \code{\link{hydromad}(sma = "gr4j", routing = "gr4jrouting")} to
#' work with models as objects (recommended).
#' @references Perrin, C., C. Michel, et al. (2003). "Improvement of a
#' parsimonious model for streamflow simulation." \emph{Journal of Hydrology}
#' 279(1-4): 275-289
#'
#' \url{http://www.cemagref.fr/webgr/Modelesgb/gr4j/fonctionnement_gr4jgb.htm}
#' @keywords models
#' @examples
#'
#' ## view default parameter ranges:
#' str(c(
#'   hydromad.getOption("gr4j"),
#'   hydromad.getOption("gr4jrouting")
#' ))
#'
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "gr4j", routing = "gr4jrouting")
#' mod0
#'
#' ## example from
#' ## http://www.cemagref.fr/webgr/Scilab/CONT_EN/HELP_HYDROGR/c_GR4J.htm
#' dat <-
#'   cbind(
#'     P = c(
#'       0, 0, 0, 0, 0, 0.04, 0.59, 0.03, 0.01, 0.16, 0.37, 8.76, 2.65,
#'       0.05, 0.02, 0.02, 0.38, 0.00, 0.02, 0.46, 4.46, 7.71, 5.71, 0.79, 1.33
#'     ),
#'     E = c(
#'       0, 0, 0, 0, 0, 0.24, 0.24, 0.24, 0.24, 0.24, 0.25, 0.25, 0.26,
#'       0.27, 0.28, 0.32, 0.33, 0.34, 0.35, 0.36, 0.36, 0.37, 0.37, 0.38, 0.38
#'     )
#'   )
#' datz <- zoo(dat, as.Date("2000-01-01") + 1:nrow(dat))
#' hydromad.options(warmup = 0)
#' modz <- hydromad(datz,
#'   sma = "gr4j", routing = "gr4jrouting",
#'   x1 = 665, x2 = 1.18, x3 = 90, x4 = 3.8, S_0 = 0.6, R_0 = 0.7
#' )
#' xyplot(predict(modz, return_state = TRUE, return_components = TRUE),
#'   strip = FALSE, strip.left = TRUE
#' )
#'
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0, x1 = 100, x2 = 20, x3 = 1, x4 = 10)
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[, 1:2], gr4j = testQ))
#'
#' ############################################
#' ## show effect of increase/decrease in each parameter
#' parRanges <- c(
#'   hydromad.getOption("gr4j")[1],
#'   hydromad.getOption("gr4jrouting")
#' )
#' parsims <- mapply(
#'   val = parRanges, nm = names(parRanges),
#'   FUN = function(val, nm) {
#'     lopar <- min(val)
#'     hipar <- max(val)
#'     names(lopar) <- names(hipar) <- nm
#'     fitted(runlist(
#'       decrease = update(mod1, newpars = lopar),
#'       increase = update(mod1, newpars = hipar)
#'     ))
#'   }, SIMPLIFY = FALSE
#' )
#'
#' xyplot.list(parsims,
#'   superpose = TRUE, layout = c(1, NA),
#'   strip = FALSE, strip.left = TRUE,
#'   main = "Simple parameter perturbation example"
#' ) +
#'   latticeExtra::layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
#'
#'
#' ############################################
#' # Example optimisation, using transformed parameters
#'
#' data(Cotter)
#' x <- Cotter[1:1000]
#'
#' # Specify gr4j model
#' mod0 <- hydromad(x, sma = "gr4j", routing = "gr4jrouting", transformed = TRUE)
#' # Use transformed parameter ranges
#' mod0 <- update(mod0, newpars = gr4j.transformpar(c(
#'   hydromad.getOption("gr4j"),
#'   hydromad.getOption("gr4jrouting")
#' )))
#' # Allow etmult to vary, because we're using temperature data instead of PET.
#' mod0 <- update(mod0, etmult = c(0.05, 1.5))
#' # Broaden a single parameter range, just as an example
#' mod0 <- update(mod0, x1 = gr4j.transformpar(list(x1 = c(100, 5000)))[["x1"]])
#'
#' mod0
#'
#' ## now try to fit the free parameters
#' set.seed(10)
#' fit1 <- fitByOptim(mod0)
#'
#' fit1
#' summary(fit1)
#' xyplot(fit1)
#'
#' # Parameters in original parameter space
#' gr4j.transformpar(coef(fit1), back = TRUE)
#' @useDynLib hydromad sma_gr4j
#' @export
gr4j.sim <-
  function(DATA,
           x1, etmult = 1, S_0 = 0.5,
           return_state = FALSE, transformed = FALSE) {
    if (transformed) x1 <- exp(x1)
    stopifnot(c("P", "E") %in% colnames(DATA))
    ## check values
    stopifnot(x1 >= 0)
    stopifnot(etmult >= 0)
    stopifnot(S_0 >= 0)
    stopifnot(S_0 <= 1)

    inAttr <- attributes(DATA[, 1])
    DATA <- as.ts(DATA)
    P <- DATA[, "P"]
    E <- DATA[, "E"] * etmult

    ## skip over missing values (maintaining the state S)
    bad <- is.na(P) | is.na(E)
    P[bad] <- 0
    E[bad] <- 0
    COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
    if (COMPILED) {
      ans <- .C(sma_gr4j,
        as.double(P),
        as.double(E),
        as.integer(NROW(DATA)),
        as.double(x1),
        as.double(S_0),
        U = double(NROW(DATA)),
        S = double(NROW(DATA)),
        ET = double(NROW(DATA)),
        NAOK = FALSE, PACKAGE = "hydromad"
      )
      U <- ans$U
      S <- ans$S
      ET <- ans$ET
    } else {
      ## implementation in R for cross-checking (slow)
      U <- S <- ET <- P
      S_prev <- S_0 * x1
      for (t in seq(1, length(P))) {
        Pn <- max(P[t] - E[t], 0)
        En <- max(E[t] - P[t], 0)
        St_x1 <- S_prev / x1
        ## production
        Ps <- 0
        ET[t] <- 0
        if (Pn > 0) {
          ## part of Pn fills the production store
          Ps <- (x1 * (1 - St_x1^2) * tanh(Pn / x1) /
            (1 + St_x1 * tanh(Pn / x1)))
        } else {
          ## actual evapo-transpiration
          ET[t] <- (S_prev * (2 - St_x1) * tanh(En / x1) /
            (1 + (1 - St_x1) * tanh(En / x1)))
        }
        S[t] <- S_prev - ET[t] + Ps
        ## percolation leakage
        perc <- S[t] * (1 - (1 + ((4 / 9) * (S[t] / x1))^4)^(-0.25))
        S[t] <- S[t] - perc
        U[t] <- perc + (Pn - Ps)
        S_prev <- S[t]
      }
    }

    attributes(U) <- inAttr
    ## re-insert missing values
    U[bad] <- NA
    ans <- U
    if (return_state) {
      attributes(S) <- attributes(ET) <- attributes(U)
      ans <- cbind(U = U, S = S, ET = ET)
    }
    return(ans)
  }


gr4j.ranges <- function() {
  list(
    x1 = c(100, 1200),
    etmult = 1
  )
}


#' @rdname gr4j
#' @useDynLib hydromad routing_gr4j
#' @export
gr4jrouting.sim <-
  function(U, x2, x3, x4, R_0 = 0, split = 0.9,
           return_components = FALSE,
           epsilon = hydromad.getOption("sim.epsilon"), transformed = FALSE) {
    if (transformed) {
      x2 <- sinh(x2)
      x3 <- exp(x3)
      x4 <- exp(x4) + 0.5
    }
    ## check values
    stopifnot(is.numeric(x2))
    stopifnot(x3 >= 0)
    stopifnot(x4 >= 0.5)
    stopifnot(R_0 >= 0)
    stopifnot(R_0 <= 1)

    inAttr <- attributes(U)
    U <- as.ts(U)

    n <- ceiling(x4)
    m <- ceiling(x4 * 2)

    ## S-curves: cumulative proportion of input with time
    n2 <- floor(x4)
    SH1 <- pmin((1:n / x4)^(5 / 2), 1)
    SH2 <- pmin(
      c(
        0 + 0.5 * (1:n2 / x4)^(5 / 2),
        1 - 0.5 * (2 - n:m / x4)^(5 / 2)
      ),
      1
    )
    SH2[1:m / x4 > 2] <- 1
    ## unit hydrographs
    UH1 <- diff(c(0, SH1))
    UH2 <- diff(c(0, SH2))

    ## skip over missing values (maintaining the state)
    bad <- is.na(U)
    U[bad] <- 0

    filter.pad0 <- function(x, f) {
      y <- x
      y[] <- filter(c(rep(0, length(f)), x),
        filter = f, sides = 1
      )[-(1:length(f))]
      y
    }

    Q9 <- filter.pad0(split * U, UH1)
    Q1 <- filter.pad0((1 - split) * U, UH2)

    COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
    if (COMPILED) {
      ans <- .C(routing_gr4j,
        as.double(Q9),
        as.double(Q1),
        as.integer(length(U)),
        as.double(x2),
        as.double(x3),
        as.double(R_0),
        Qr = double(length(U)),
        Qd = double(length(U)),
        R = double(length(U)),
        NAOK = FALSE, PACKAGE = "hydromad"
      )
      Qr <- ans$Qr
      Qd <- ans$Qd
      R <- ans$R
    } else {
      ## implementation in R for cross-checking (slow)
      Qd <- Qr <- R <- U
      R_prev <- R_0 * x3
      for (t in seq(1, length(U))) {
        Rt_x3 <- R_prev / x3
        ## groundwater exchange term
        F <- x2 * Rt_x3^(7 / 2)
        ## reservoir level
        R[t] <- max(0, R_prev + Q9[t] + F)
        ## outflow of reservoir
        Qr[t] <- R[t] * (1 - (1 + (R[t] / x3)^4)^(-0.25))
        R[t] <- R[t] - Qr[t]
        ## other store
        Qd[t] <- max(0, Q1[t] + F)
        R_prev <- R[t]
      }
    }

    ## zap simulated values smaller than epsilon
    Qr[Qr < epsilon] <- 0
    Qd[Qd < epsilon] <- 0

    attributes(Qr) <- attributes(Qd) <- attributes(R) <- inAttr

    ## re-insert missing values
    Qr[bad] <- NA
    Qd[bad] <- NA

    if (return_components) {
      return(cbind(Xr = Qr, Xd = Qd, R = R))
    } else {
      return(Qr + Qd)
    }
  }


gr4jrouting.ranges <- function() {
  list(
    x2 = c(-5, 3),
    x3 = c(20, 300),
    x4 = c(1.1, 2.9)
  )
}



#' Transform GR4J parameters
#'
#' Apply or reverse transformation of GR4J parameters
#'
#'
#' @param pars named vector or list of parameters, e.g. as provided by
#' \code{coef.hydromad}.
#' @param back Whether to transform or untransform (reverse) the parameters.
#' @return Named list of transformed/untransformed parameters, depending on
#' value of \code{back}.
#' @author Joseph Guillaume
#' @seealso \code{\link{gr4j}}
#' @keywords models
#' @examples
#'
#' gr4j.transformpar(c(hydromad.getOption("gr4j"), hydromad.getOption("gr4jrouting")))
#' gr4j.transformpar(c(x1 = 150, x2 = 2, x3 = 50, x4 = 2), back = FALSE)
#' @export
gr4j.transformpar <- function(pars, back = FALSE) {
  pars <- modifyList(list(x1 = NA, x2 = NA, x3 = NA, x4 = NA), as.list(pars))
  newpars <- pars
  if (back) {
    newpars[["x1"]] <- exp(pars[["x1"]])
    newpars[["x2"]] <- sinh(pars[["x2"]])
    newpars[["x3"]] <- exp(pars[["x3"]])
    newpars[["x4"]] <- exp(pars[["x4"]]) + 0.5
  } else {
    newpars[["x1"]] <- log(pars[["x1"]])
    newpars[["x2"]] <- asinh(pars[["x2"]])
    newpars[["x3"]] <- log(pars[["x3"]])
    newpars[["x4"]] <- log(pars[["x4"]] - 0.5)
  }
  newpars
}
