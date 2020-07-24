## poweroid geometry
## (c) Felix Andrews <felix@nfrac.org> 2009-06-04

## poweroid(V, alpha, beta) --> V, A, H, r
## poweroid(A, alpha, beta) --> V, A, H, r
## poweroid(H, alpha, beta) --> V, A, H, r
## poweroid(r, alpha, beta) --> V, A, H, r
## poweroid(V, A, alpha) --> beta
## poweroid(V, A, beta) --> alpha
## poweroid(V, r, alpha) --> beta
## poweroid(V, r, beta) --> alpha
## poweroid(V[:], A[:]) --> alpha, beta (fit)



#' Poweroid geometry (cones, paraboloids, etc).
#'
#' Poweroid geometry (cones, paraboloids, etc).
#'
#' A poweroid is a generalisation of a cone, paraboloid (3D parabola), etc.
#'
#' A poweroid is defined by the surface of rotation (around the height axis) of
#' the curve: \deqn{H = \alpha * r ^ \beta}{H = alpha * r ^ beta}
#'
#' where H is height, r is radius, and alpha and beta are parameters.
#'
#' Then, as special cases, \describe{ \item{list("beta = 1")}{ corresponds to a
#' cone. } \item{list("beta = 2")}{ corresponds to a paraboloid. }
#' \item{list("beta = 4")}{ is sometimes called a quartoid. } \item{list("beta
#' < 1")}{ has a funnel shape. } }
#'
#' The \code{poweroid} function can be used in three different ways:
#'
#' If \code{alpha} and \code{beta} are specified, any one of \code{V},
#' \code{H}, \code{A} or \code{r} can be used to calculate all the others.
#'
#' If \code{beta} is specified, any two of \code{V}, \code{H}, or \code{A / r}
#' can be used to derive \code{alpha}.
#'
#' Otherwise, any two of \code{V}, \code{H}, or \code{A / r} can be used to
#' estimate \code{alpha} and \code{beta}. In this case, multiple values must be
#' given for each measurement, and the parameters are fitted to these using
#' \code{\link{optimize}}. (Specifically, the beta parameter is fitted, while
#' the alpha parameter is derived from beta on each iteration.)
#'
#' The examples below demonstrate the different shapes and depth distributions
#' that are possible.
#'
#' @importFrom stats weighted.mean nlminb optimize
#'
#' @param alpha scale parameter.
#' @param beta shape parameter.
#' @param V volume.
#' @param H height.
#' @param A area (of circular base).
#' @param r radius (of circular base).
#' @param \dots passed to \code{\link{optimize}} when fitting parameters.  The
#' following arguments apply only to the case of fitting parameters.
#' @param rel.error if \code{TRUE}, consider errors as proportion of given
#' observed values; otherwise consider the total error in original units. This
#' affects the objective function as well as estimation of alpha.
#' @param polish to apply \code{\link{nlminb}} when fitting parameters.
#' @param details to return calibration details when fitting parameters.
#' @return If \code{alpha} and \code{beta} are specified, a \code{data.frame}
#' with columns \code{V}, \code{H}, \code{A} and \code{r}.
#'
#' Otherwise, a list with components \code{alpha} and \code{beta}.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{swimp}}, \code{\link{convertFlow}}
#' @references \url{http://mathworld.wolfram.com/Paraboloid.html}
#'
#' \url{http://mathworld.wolfram.com/Poweroid.html}
#' @keywords math optimize
#' @examples
#'
#' ## USE CASE 1: specified parameters
#'
#' ## A cone with height H and slope alpha
#' poweroid(H = 1:10, alpha = 2, beta = 1)
#'
#' ## Take four different shapes: beta = 4 / 2 / 1 / 0.5
#' ## Calculate geometry, given Volume values from 0 to 100.
#' shapedat <- list(
#'   `beta=4 (quartoid)` = poweroid(V = (0:1000) / 10, beta = 4, alpha = 0.0078125),
#'   `beta=2 (paraboloid)` = poweroid(V = (0:1000) / 10, beta = 2, alpha = 0.125),
#'   `beta=1 (cone)` = poweroid(V = (0:1000) / 10, beta = 1, alpha = 0.5),
#'   `beta=1/2 (funnel)` = poweroid(V = (0:1000) / 10, beta = 0.5, alpha = 1)
#' )
#'
#' ## Plot cross sections, all with V = 100 and intersecting at r = 4:
#'
#' shapedf <- do.call(make.groups, shapedat)
#'
#' xyplot(c(H, H) ~ c(-r, r), shapedf,
#'   groups = rep(which, 2), type = "a", lwd = 2,
#'   xlab = "radius", ylab = "height",
#'   auto.key = list(lines = TRUE, points = FALSE)
#' )
#'
#' ## Compare general form of depth distributions:
#'
#' depthdat <- lapply(shapedat, transform,
#'   H = (1 - H / max(H)), A = 100 * A / max(A)
#' )
#' depthdf <- do.call(make.groups, depthdat)
#'
#' xyplot(H ~ A | which, depthdf,
#'   as.table = TRUE,
#'   xlab = "cumulative area (percent)", ylab = "depth / max(depth)",
#'   panel = function(x, y, ...) {
#'     panel.grid(h = -1, v = -1)
#'     panel.polygon(c(x, 0), c(y, 0),
#'       col = grey(0.75),
#'       border = "transparent"
#'     )
#'     panel.xyplot(x, y, ...)
#'   },
#'   type = "l", par.settings = simpleTheme(lwd = 2)
#' )
#'
#' ## Compare effect of different alpha, with a fixed Volume V = 100
#' ## beta is constant at 0.8
#'
#' scaledat <- list(
#'   `alpha = 2` = poweroid(V = (0:1000) / 10, beta = 0.8, alpha = 2),
#'   `alpha = 1` = poweroid(V = (0:1000) / 10, beta = 0.8, alpha = 1),
#'   `alpha = 1/2` = poweroid(V = (0:1000) / 10, beta = 0.8, alpha = 0.5),
#'   `alpha = 1/4` = poweroid(V = (0:1000) / 10, beta = 0.8, alpha = 0.25)
#' )
#'
#' scaledf <- do.call(make.groups, scaledat)
#'
#' xyplot(c(H, H) ~ c(-r, r), scaledf,
#'   groups = rep(which, 2), type = "a", lwd = 2,
#'   xlab = "radius", ylab = "height",
#'   auto.key = list(lines = TRUE, points = FALSE)
#' )
#'
#' depthdat <- lapply(scaledat, transform, H = max(H) - H)
#' depthdf <- do.call(make.groups, depthdat)
#'
#' xyplot(H ~ A,
#'   groups = which, depthdf,
#'   xlab = "cumulative area", ylab = "depth",
#'   panel = function(...) {
#'     panel.superpose(..., panel.groups = function(x, y, ...) {
#'       panel.polygon(c(x, 0), c(y, 0), col = grey(0.25), alpha = 0.25)
#'     })
#'     panel.xyplot(...)
#'   },
#'   type = "l", par.settings = simpleTheme(lwd = 2),
#'   auto.key = list(lines = TRUE, points = FALSE)
#' )
#'
#' ## USE CASE 2: Derive alpha given beta.
#'
#' poweroid(V = 1000, A = 20, beta = 2)
#' poweroid(V = 1000, H = 5, beta = 1.5)
#' ## derive the alpha values used at the top of this section
#' poweroid(r = 4, H = 2, beta = c(4, 2, 1, 0.5))
#'
#' ## USE CASE 3: fit alpha and beta to given data
#'
#' poweroid(V = c(100, 1000, 10000), A = c(5, 52, 127))
#' # --> alpha = 4.83, beta = 1.9
#' ## simulate with these parameters to see how well they fit A:
#' poweroid(V = c(100, 1000, 10000), alpha = 4.83, beta = 1.9)
#' # --> A:  11.95  38.91 126.70
#'
#' ## try fitting again, based on relative error this time:
#' poweroid(V = c(100, 1000, 10000), A = c(5, 52, 127), rel.error = TRUE)
#' # --> alpha = 26.6, beta = 1.14
#' ## simulate with these parameters to see how well they fit A:
#' poweroid(V = c(100, 1000, 10000), alpha = 26.6, beta = 1.14)
#' # --> A:   6.72  29.11 126.20
#' @export
poweroid <-
  function(alpha = NULL, beta = NULL,
           V = NULL, H = NULL,
           A = if (!missing(r)) pi * r^2,
           r = if (!missing(A)) sqrt(A / pi),
           ...,
           rel.error = FALSE,
           polish = FALSE, details = FALSE) {
    force(list(A, r))
    objective <-
      if (rel.error) {
        ~ mean(abs((obs - mod) / obs))
      } else {
        ~ mean(abs((obs - mod)))
      }
    if (!is.null(alpha) && !is.null(beta)) {
      ## geometry specified parametrically
      if (!is.null(V)) {
        A <- ((2 / beta + 1) * V / (alpha * pi^(-beta / 2)))^(2 / (2 + beta))
      } else if (!is.null(H)) {
        ## NOTE: this form may blow up with beta < ~0.01
        r <- (H / alpha)^(1 / beta)
        A <- pi * r^2
      }
      if (is.null(A)) {
        stop("give one of V, H, A or r")
      }
      ## now we have A, derive the rest
      r <- sqrt(A / pi)
      V <- (alpha * pi^(-beta / 2) / (2 / beta + 1)) * A^(1 + beta / 2)
      H <- alpha * r^beta
      return(data.frame(V = V, A = A, H = H, r = r))
    }
    if (!is.null(beta)) {
      ## beta (shape type) specified, derive alpha exactly
      if (!is.null(H) && !is.null(r)) {
        alpha <- H / (r^beta)
      } else
      if (!is.null(V) && !is.null(r)) {
        alpha <- (2 / beta + 1) * V / (pi * r^(2 + beta))
      } else
      if (!is.null(V) && !is.null(H)) {
        alpha <- (pi / ((2 / beta + 1) * V))^(beta / 2) * H^(1 + beta / 2)
        ## should avoid H ^ 2/beta
        ## alpha <- ((pi * H ^ (2/beta + 1)) / ((2 / beta + 1) * V)) ^ (beta/2)
      } else {
        stop("give any two of V, H, A/r to derive alpha")
      }
      return(list(alpha = alpha, beta = beta))
    }

    ## shape not specified; fit it to given data

    estAlpha <- function(betaEst) {
      ## find alpha values that fit data exactly, given this beta
      ## then take a (weighted) average
      alphaVec <- poweroid(V = V, H = H, A = A, beta = betaEst)$alpha
      if (rel.error) {
        alpha <- mean(alphaVec)
      } else {
        w <- if (!is.null(A)) A else H
        alpha <- weighted.mean(alphaVec, w = w)
      }
      alpha
    }
    objfn <- function(x, withAlpha = FALSE) {
      if (withAlpha) {
        alpha <- x[1]
        beta <- x[2]
      } else {
        beta <- x[1]
        alpha <- estAlpha(beta)
      }
      if (!is.null(H) && !is.null(A)) {
        H.hat <- alpha * r^beta
        foo <- list(obs = H, mod = H.hat)
      } else
      if (!is.null(V) && !is.null(A)) {
        A.hat <- ((2 / beta + 1) * V / (alpha * pi^(-beta / 2)))^(2 / (2 + beta))
        foo <- list(obs = A, mod = A.hat)
      } else
      if (!is.null(V) && !is.null(H)) {
        H.hat <- ((2 / beta + 1) * V / pi)^(beta / (2 + beta)) * alpha^(2 / (2 + beta))
        foo <- list(obs = H, mod = H.hat)
      } else {
        stop("need more information")
      }
      objtmp <- objective
      if (inherits(objective, "formula")) {
        objtmp <- objective[[2]]
      }
      eval(objtmp, foo)
    }
    ## run optimize with different starting points for beta
    results <- lapply(
      c(0.5, 1, 1.5, 2, 2.5, 3.5),
      function(beta) {
        optimize(objfn, interval = pmax(0.01, beta + c(-0.5, 0.5)), ...)
      }
    )
    if (details) {
      return(results)
    }
    ## extract overall best result (value of beta)
    whichBest <- which.min(sapply(results, function(x) x$objective))
    beta <- results[[whichBest]]$minimum
    alpha <- estAlpha(beta)
    ## now polish result
    ## this does not seem to change the result at all:
    if (polish) {
      result <- nlminb(c(alpha = alpha, beta = beta), objfn,
        lower = c(0, 0.01),
        withAlpha = TRUE
      )
      if (result$convergence != 0) {
        warning(result$message)
      }
      alpha <- result$par[1]
      beta <- result$par[2]
    }
    return(list(alpha = alpha, beta = beta))
    #    result <- nlm(objfn, c(alpha = alpha, beta = beta),
    #                     withAlpha = TRUE, ...)
    #    if (result$code == 4) {
    #        warning("iteration limit exceeded")
    #    } else if (result$code > 3) {
    #        warning("code ", result$code, ", see ?nlm")
    #    }
    #    return(list(alpha = result$estimate[1],
    #                beta = result$estimate[2]))
  }
