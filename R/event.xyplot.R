## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


#' Scatterplots with variables aggregated in event windows
#'
#' Scatterplots with variables aggregated in event windows.
#'
#' @importFrom stats terms model.frame model.response as.formula
#' @importFrom utils stack
#'
#' @name event.xyplot
#' @aliases event.xyplot.formula
#' @param x an object for the generic method; in the \code{formula} method this
#' is a \code{\link{formula}} as described in the Details section.  All terms
#' in the formula should be aggregated in the given \code{events}; for this
#' purpose a special function \code{e()} is defined as a shorthand for
#' \code{\link{eventapply}}, with arguments \code{e(X, FUN = sum, ...)}.  The
#' formula may specify a response (the left hand side of formula), which is
#' treated as a single term, although it can be a matrix with multiple columns.
#' Alternatively, the response matrix can be given in \code{response}.  The
#' covariate terms are also allowed to be matrix-like with the same number of
#' columns as the response, in which case each column is plotted against the
#' corresponsing column of the response.
#' @param data data source giving variables used in the formula. Typically a
#' \code{data.frame}. Passed to \code{\link{model.frame}}.
#' @param events event sequence produced by \code{\link{eventseq}}, or a vector
#' defining continguous groups on the specified variables.
#' @param \dots further arguments passed to \code{\link{xyplot}} (actually
#' \code{\link{xyplot.list}}).
#' @param response a response vector or matrix corresponding to the variables
#' in the formula. If this is given, any response in the formula is ignored.
#' @param eFUN a function to aggregate \code{response}, if it is given, in the
#' given \code{events}. If a reponse appears in the formula it should be in a
#' \code{e()} term, and \code{eFUN} is not used.
#' @param as.table,layout,auto.key,xlab,ylab passed to \code{\link{xyplot}}.
#' @return this function returns a trellis object which can be \code{plot}ted.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{xyplot.list}}, \code{\link{xyplot}},
#' \code{\link{panel.xyplot}}, \code{\link{eventapply}}
#' @keywords hplot ts
#' @examples
#'
#' data(Canning)
#' ev <- eventseq(Canning$P,
#'   thresh = 20, inthresh = 1, indur = 3,
#'   continue = FALSE
#' )
#'
#' event.xyplot(e(Q, sum) / e(P, sum) ~ e(P, sum) + e(P, max) + e(lag(Q, -1), first),
#'   data = Canning, events = ev,
#'   scales = list(y = list(log = TRUE)),
#'   yscale.components = yscale.components.fractions,
#'   ylab = "event runoff ratio (Q/P)", layout = c(3, 1),
#'   xlab = c(
#'     "total rain (mm)", "max rain (mm/day)",
#'     "antecedent flow (mm/day)"
#'   )
#' )
#'
#' ## multiple response variables
#' event.xyplot(cbind(e(Q, quantile, 0.9), e(Q, quantile, 0.99)) ~
#' e(P, max) + e(P, mean) + e(P, sum),
#' data = Canning, events = ev,
#' scales = list(y = list(log = TRUE)),
#' yscale.components = yscale.components.log10.3
#' )
#' @export event.xyplot
event.xyplot <- function(x, ...) {
  UseMethod("event.xyplot")
}

#' @rdname event.xyplot
#' @export
event.xyplot.formula <-
  function(x, data = list(), events, ...,
           response = NULL, eFUN = sum,
           as.table = TRUE, layout = NULL,
           auto.key = NCOL(response) > 1,
           xlab = "event-aggregated covariate values",
           ylab = "event-aggregated response") {
    x <- as.formula(x)
    ## for zoo objects this splits up the columns into items:
    data <- as.list(data)
    ## define a special function to aggregate events, for use in formula
    e <- function(X, FUN = sum, ...) eventapply(X, events, FUN = FUN, ...)
    ## make it known to the formula:
    environment(x) <- environment()
    ## another convenience function for formula:
    first <- function(x) head(x, 1)
    ## extract covariates from formula
    tt <- terms(x, data = data, keep.order = TRUE)
    mf <- model.frame(tt, data, na.action = na.pass) ## na.omit stuffs up index
    if (is.null(response)) {
      response <- model.response(mf)
      if (is.null(response)) {
        stop("no response in formula, and 'response' argument missing")
      }
      mf <- mf[, -1, drop = FALSE]
    } else {
      response <- e(response, eFUN)
    }
    ystack <- stack(as.data.frame(response))
    foo <-
      xyplot.list(mf,
        data = ystack,
        FUN = function(VAR, ...) {
          xyplot(values ~ rep(coredata(VAR), length = NCOL(response) * NROW(response)),
            groups = ind, ...
          )
        },
        default.scales = list(x = list(relation = "free")),
        ...,
        as.table = as.table, layout = layout,
        auto.key = auto.key,
        xlab = xlab, ylab = ylab,
        y.same = NA
      ) ## allow specified 'ylim'
    foo$call <- sys.call(sys.parent())
    foo
  }

#' @import utils
utils::globalVariables(c("ind"))
