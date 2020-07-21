## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##




#' Compare calibrations with different transfer function (ARMA) orders for
#' routing.
#'
#' Compare calibrations with different transfer function (ARMA) orders for
#' routing.
#'
#'
#' @aliases tryModelOrders summary.tryModelOrders
#' @param expr an expression to calibrate a hydromad model.
#' @param n the set of values of \var{n} to try (order of the auto-regressive
#' component).
#' @param m the set of values of \var{m} to try (order of the moving-average
#' component).
#' @param delay the set of delays to try.
#' @param verbose show detailed tracing output.
#' @return a list of model objects, of class \code{hydromad.\link{runlist}}.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{armax}}, \code{\link{armax.sriv.fit}}
#' @references P C Young?
#' @keywords optimize
#'
#'
#' @export
tryModelOrders <-
  function(expr, n = 0:3, m = 0:2,
           delay = hydromad.getOption("delay"),
           verbose = hydromad.getOption("trace")) {
    expr <- substitute(expr)
    allOrders <- expand.grid(m = m, n = n, delay = delay)
    allOrders <- subset(allOrders, n >= m)
    results <- list()
    beSilent <- hydromad.getOption("quiet")
    ## restore previous setting when finished
    oopt <- hydromad.options("order", "delay", "quiet")
    on.exit(hydromad.options(oopt))
    hydromad.options(quiet = !verbose)
    for (i in 1:NROW(allOrders)) {
      orderSpec <- allOrders[i, ]
      order <- c(n = orderSpec$n, m = orderSpec$m)
      d <- orderSpec$delay
      hydromad.options(order = order, delay = d)
      if (any(!is.na(delay))) {
        hydromad.options(delay = d)
        nm <- paste("(n=", order[1], ", m=", order[2],
          ", d=", d, ")",
          sep = ""
        )
      } else {
        nm <- paste("(n=", order[1], ", m=", order[2], ")", sep = "")
      }
      if (!isTRUE(beSilent)) {
        message("order: ", nm, "... ", appendLF = FALSE)
      }
      ## execute the expression and store the result
      mod <- eval.parent(expr)
      results[[nm]] <- mod
      if (!isTRUE(beSilent)) {
        if (!isFullySpecified(mod)) {
          message("model not fully specified!")
        } else if (!isValidModel(mod)) {
          message(toString(mod, width = 60))
        } else {
          modsumm <- summary(mod,
            stats = c("ARPE", "r.squared"),
            with.hydrostats = FALSE
          )
          ARPE <- modsumm$ARPE
          R2 <- modsumm$r.squared
          message(
            " ARPE = ", round(ARPE, 3),
            ", R^2 = ", round(R2, 3)
          )
        }
      }
    }
    ans <- as.runlist(results)
    class(ans) <- c("tryModelOrders", class(ans))
    ans
  }


#' @export
summary.tryModelOrders <-
  function(object,
           stats = c("ARPE", "r.squared", "r.sq.log"),
           ...) {
    class(object) <- setdiff(class(object), "tryModelOrders")
    summary(object, stats = stats, ...)
  }
