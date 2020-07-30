## hydromad: Hydrological Modelling and Analysis of Data

#' Simple degree day factor snow model
#'
#' Simple degree day factor snow model coupled with IHACRES CMD soil moisture
#' model.
#'
#' SWE snow water equivalent
#'
#' ISWE water equivalent of ice in the snowpack
#'
#' LSWE liquid water retained in the snowpack
#'
#'
#' @name snow
#' @aliases snow.sim
#' @param DATA a \code{\link{ts}}-like object with named columns: \describe{
#' \item{list("P")}{ time series of areal rainfall depths, usually in mm. }
#' \item{list("E")}{ time series of potential evapo-transpiration, or more
#' typically, temperature as an indicator of this. } }
#' @param Tmax temperature threshold for rain, all rain is liquid above this
#' threshold.
#' @param Tmin temperature threshold for rain, all rain is snow below this
#' threshold.
#' @param Tmelt temperature threshold for snowmelt and freezing in the
#' snowpack.
#' @param kd degree day factor for snowmelt.
#' @param kf degree day factor for freezing.
#' @param rcap retention parameter for liquid water capacity of snowpack.
#' @param cr correction factor for rainfall.
#' @param cs correction factor for snowfall.
#' @param LSWE_0,ISWE_0 initial values of state variables.
#' @param \dots parameters for the \link{IHACRES.CMD.model}.
#' @param return_state to return state variables as well as the effective
#' rainfall.
#' @return \code{snow.sim} returns the modelled time series of effective
#' rainfall, or if \code{return_state = TRUE}, a multi-variate time series with
#' named columns \code{U} (effective rainfall), \code{SWE} (snow water
#' equivalent) and \code{TF}, as well as the CMD state variables.
#' @author Coded in R by Jarkko Koskela @@tkk.fi 2010-02-26.
#'
#' Converted to C by Felix Andrews \email{felix@@nfrac.org}.
#' @seealso \code{\link{hydromad}(sma = "snow")} to work with models as objects
#' (recommended).
#' @references Kokkonen T., Jakeman A.J, Koivusalo.H, Norton.J.: COMPUTATIONAL
#' METHODS FOR WATER RESOURCE ASSESSMENTS: AN EXERCISE KIT Educational Series
#' on Modelling and Software iEMSs International Modelling and Software Society
#' Available through www.iemss.org
#' @keywords models
#' @examples
#'
#' ## view default parameter ranges:
#' str(hydromad.options("snow"))
#'
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "snow", routing = "expuh")
#' mod0
#'
#' ## simulate with some arbitrary parameter values
#' mod1 <- update(mod0,
#'   Tmax = 15, Tmin = 5, cr = 1, cs = 1,
#'   kd = 3, kf = 1, rcap = 0.5,
#'   d = 200, f = 0.5, e = 0.1, tau_s = 10
#' )
#'
#' ## plot results with state variables
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(cbind(HydroTestData[, 1:2], snow = testQ))
#'
#' ## show effect of increase/decrease in each parameter
#' parlist <- list(
#'   Tmax = c(10, 20), Tmin = c(0, 10),
#'   cr = c(0.5, 2), cs = c(0.5, 2),
#'   kd = c(2, 5), kf = c(0, 2), rcap = c(0, 1)
#' )
#' parsims <- mapply(
#'   val = parlist, nm = names(parlist),
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
#' @useDynLib hydromad sma_snow
#' @export
snow.sim <-
  function(DATA, Tmax, Tmin, kd, kf, rcap, Tmelt = Tmin,
           cr = 1, cs = 1, LSWE_0 = 0, ISWE_0 = 0,
           ..., return_state = FALSE) {
    stopifnot(c("P", "E") %in% colnames(DATA))
    ## check values
    stopifnot(0 <= kd)
    stopifnot(0 <= kf)
    stopifnot(0 <= cr)
    stopifnot(0 <= cs)
    stopifnot(0 <= rcap)
    Tmin <- min(Tmax, Tmin)

    inAttr <- attributes(DATA[, 1])
    DATA <- as.ts(DATA)
    P <- DATA[, "P"]
    E <- DATA[, "E"]

    ## rainfall or snowfall
    fr <- (E - Tmin) / (Tmax - Tmin)
    fr <- pmax(pmin(fr, 1), 0)
    Prain <- fr * cr * P
    Psnow <- (1 - fr) * cs * P

    COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
    if (COMPILED) {
      ans <- .C(sma_snow,
        as.double(Prain),
        as.double(Psnow),
        as.double(E),
        as.integer(NROW(DATA)),
        as.double(kd),
        as.double(kf),
        as.double(rcap),
        as.double(Tmelt),
        as.double(LSWE_0),
        as.double(ISWE_0),
        U = double(NROW(DATA)),
        SWE = double(NROW(DATA)),
        NAOK = FALSE, PACKAGE = "hydromad"
      )
      Sdischarge <- ans$U
      SWE <- ans$SWE
    } else {
      ## implementation in R for cross-checking (slow)
      ## time loop
      SWE <- Sdischarge <- P * 0
      LSWEprev <- LSWE_0
      ISWEprev <- ISWE_0
      for (t in seq(1, length(P))) {
        ## Melt (degree day model)
        melt <- min(max(kd * (E[t] - Tmelt), 0), ISWEprev)
        ## Freezing (degree day model)
        freeze <- min(max(kf * (Tmelt - E[t]), 0), LSWEprev)
        ## Mass balance for the snowpack
        ##
        ## Ice in the snowpack
        ISWE <- ISWEprev + Psnow[t] + freeze - melt
        ## Water in the snowpack
        LSWE <- min(rcap * ISWE, LSWEprev + Prain[t] + melt - freeze)
        ## Rain/melt is snowmelt discharge when there is snow on the ground,
        ## and rainfall in snow-free periods.
        Sdischarge[t] <- max(Prain[t] + melt - freeze - (rcap * ISWE - LSWEprev), 0)
        SWE[t] <- LSWE + ISWE
        ISWEprev <- ISWE
        LSWEprev <- LSWE
      }
    }

    DATA[, "P"] <- Sdischarge

    ## IHACRES CMD-module
    U <- cmd.sim(DATA, ..., return_state = return_state)

    if (return_state) {
      cmdState <- U
      attributes(cmdState) <-
        modifyList(inAttr, list(
          dim = dim(cmdState),
          dimnames = dimnames(cmdState)
        ))
      attributes(SWE) <- attributes(Sdischarge) <- inAttr
      ans <- cbind(cmdState, SWE = SWE, TF = Sdischarge)
      #        colnames(ans)[1:NCOL(cmdState)] <- colnames(U) ## TODO
    } else {
      attributes(U) <- inAttr
      ans <- U
    }
    return(ans)
  }

snow.ranges <- function() {
  list(
    Tmax = c(0, 2),
    Tmin = c(-1, 1),
    cr = c(0.8, 2),
    cs = c(0.8, 2),
    kd = c(2, 5),
    kf = c(0, 2),
    rcap = c(0, 1),
    f = c(0.01, 3),
    e = c(0.01, 1.5),
    d = 200
  )
}
