## hydromad: Hydrological Modelling and Analysis of Data
## SimHydrouting
## Coded based on diagram and description in:
# Chiew et al 2009 WATER RESOURCES RESEARCH, VOL. 45, W10414, doi:10.1029/2008WR007338, 2009

#' Routing based on Muskingum
#' 
#' Coded based on diagram and description in:
#' Chiew et al 2009 WATER RESOURCES RESEARCH, 
#' VOL. 45, W10414, doi:10.1029/2008WR007338, 2009


#' @name simhydrouting
#' @aliases simhydrouting.sim simhydrouting.ranges
#'
#' @param U effective rainfall series from a SMA object
#' @param DELAY delay parameter in the Muskingum routing (Days)
#' @param X storage weight parameter in Muskingum routing (-)
#' @param return_components included here to align with hydromad structure
#' @author Willem Vervoort \email{willemvervoort@@gmail.com} and Joseph Guillaume
#' \email{josephguillaume@@gmail.com}
#' @seealso \code{\link{hydromad}(sma = "simhyd", routing = "simhydrouting")} to
#' work with models as objects (recommended).
#' @references Chiew, F. H. S., Teng, J., Vaze, J., Post, D. A., Perraud, J. M.,
#' Kirono, D. G. C., and Viney, N. R. (2009), 
#' Estimating climate change impact on runoff across southeast Australia: 
#' Method, results, and implications of the modeling method, 
#' \emph{Water Resour. Res.}, 45, W10414, \url{doi:10.1029/2008WR007338}.
#'
#' @keywords models
#' @examples
#'
#' ## view default parameter ranges:
#' str(hydromad.getOption("simhydrouting"))
#' 
#'
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "simhyd", routing = "simhydrouting")
#' mod0
#'


#' @export
simhydrouting.sim <- function(U, DELAY = 1, X_m = 0.2,
                              epsilon = hydromad.getOption("sim.epsilon"),
                              return_components = FALSE) {
  X <- rep(0, length(U))
  inAttr <- attributes(U)
  U <- as.ts(U)
  bad <- is.na(U)
  U[bad] <- 0
  if (2 * DELAY * X_m < 1 & 2 * DELAY * (1 - X_m) > 1) {
    # Muskingum components
    C0 <- (-DELAY * X_m + 0.5) / (DELAY * (1 - X_m) + 0.5)
    # print(C0)
    C1 <- (DELAY * X_m + 0.5) / (DELAY * (1 - X_m) + 0.5)
    # print(C1)
    C2 <- (DELAY * (1 - X_m) - 0.5) / (DELAY * (1 - X_m) + 0.5)
    # print(C2)
  } else {
    C0 <- 0
    C1 <- 1
    C2 <- 0
    # print("model parameters adjusted")
  }
  
  if (C0 + C1 + C2 != 1) {
    C0 <- 0
    C1 <- 1
    C2 <- 0
    # print("model parameters adjusted again")
  }
  # print(C0+C1+C2)
  # if (round(C0+C1+C2)!=1)  C0 <- 0; C1 <- 1; C2 <- 0
  
  X[1] <- U[1]
  for (t in 1:(length(U) - 1)) {
    X[t + 1] <- C0 * U[t + 1] + C1 * U[t] + C2 * X[t]
    # print(X[t+1])
  }
  X[abs(X) < epsilon] <- 0
  X[bad] <- NA
  attributes(X) <- inAttr
  X
}

#' @export
simhydrouting.ranges <- function() {
  list(
    DELAY = c(0.1, 5),
    X_m = c(0.01, 0.5)
  )
}

