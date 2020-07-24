#' Title placeholder
#'
#' Description placeholder.
#'
#' @name eval_fdc
#'
#' @importFrom stats qnorm pnorm
#'
#' @aliases fdc.sample, fdc.allpoints
#' @param Q Placeholder
#'
#' @export
fdc.sample <- function(Q) {
  if (is.zoo(Q)) Q <- coredata(Q)
  pp <- qnorm(ppoints(NROW(Q)))
  space <- diff(pp[1:2])
  eval.prob <- pnorm(seq(min(pp), max(pp), by = space))
  fdc.q <- quantile(Q, 1 - eval.prob)
  invisible(fdc.q)
}

#' @export
fdc.allpoints <- function(Q) quantile(Q, 1 - ppoints(NROW(Q)))


# hmadstat("r.squared")(Q=fdc.sample(Murrindindi$Q),
#                      X=fdc.sample(Murrindindi$Q+1))
