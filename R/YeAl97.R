#' Model performance statistics from Ye et al. 1997
#'
#' Performance statistics from calibrations of three model structures in three
#' catchments in two data periods, assessed using absolute mean deviation, bias
#' and efficiency.
#'
#'
#' @name YeAl97
#' @docType data
#' @format A data frame with 36 observations on the following 7 variables.
#' \describe{ \item{Catchment}{a character vector, either Canning,
#' Salmon or Stones} \item{calib.period}{a character vector, either
#' "First 5Y" or "Second 5Y"} \item{sim.period}{a character vector,
#' either "First 5Y" or "Second 5Y"} \item{Model.str}{a factor with
#' levels \code{GSFB} \code{IHACRES} \code{LASCAM}} \item{A}{a numeric
#' vector, absolute mean deviation} \item{B}{a numeric vector, bias}
#' \item{E}{a numeric vector, efficiency} }
#' @references Used by \code{\link{paretoTimeAnalysis}}.
#' @source Ye, W., B.C. Bates, N.R. Viney, M. Sivapalan and A.J. Jakeman
#' (1997). Performance of conceptual rainfall-runoff models in low-yielding
#' ephemeral catchments. \emph{Water Resour. Res.} 33(1): 153-166 DOI:
#' \url{http://dx.doi.org/10.1029/96wr02840}.
#' @keywords datasets
#' @examples
#'
#' data(YeAl97)
#' YeAl97
"YeAl97"
