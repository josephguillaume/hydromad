## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##



#' Work with a set of model runs
#'
#' A \code{runlist} object is simply a list of model objects.
#'
#' Note that the \code{coef} method just calls \code{summary(..., FUN = coef)}.
#'
#' @importFrom parallel parLapply
#'
#' @aliases runlist as.runlist [.runlist c.runlist coef.runlist summary.runlist
#' print.summary.runlist residuals.runlist fitted.runlist
#' @param \dots for \code{runlist}, a named list of model objects, specified
#' directly as in \code{\link{list}}.
#'
#' In other cases, arguments are passed on to the generic functions.
#' @param x a simple \code{list}.
#' @param object a \code{runlist}: a list of fitted model objects.
#' @param FUN function returning one or more named numeric values. Any returned
#' values other than single numeric values will be ignored.
#' @param items if given, this is used to extract elements of the result from
#' \code{FUN}; otherwise, all single numeric elements are extracted.
#' @param recursive Placeholder
#' @return \code{runlist} and \code{as.runlist} return a list of class
#' \code{runlist} and also (firstly) \var{class}\code{.runlist}, where
#' \var{class} is the first class of the first element of the list.
#'
#' \code{summary} and \code{coef} return a data frame, with rows for each
#' element of \code{object} and columns for each named item returned by
#' \code{FUN}. Any missing items will filled in with \code{NA}.
#' @author Felix Andrews \email{felix@@nfrac.org}
#' @seealso \code{\link{xyplot.runlist}}, \code{coef.hydromad},
#' \code{\link{summary.hydromad}}
#' @keywords utilities
#'
#'
#'
#' @export
runlist <- function(...) {
  object <- list(...)
  if (is.null(names(object))) {
    names(object) <- rep("", length(object))
  }
  unnamed <- names(object) == ""
  if (any(unnamed)) {
    dnames <- as.list(substitute(list(...)))[-1]
    for (i in seq_along(object)) {
      if (names(object)[i] == "") {
        names(object)[i] <- toString(deparse(dnames[[i]]), width = 12)
      }
    }
  }
  names(object) <- make.unique(names(object))
  if (length(object) > 0) {
    class(object) <- paste(class(object[[1]]), "runlist", sep = ".")
  }
  class(object) <- unique(c(class(object), "runlist", "list"))
  object
}


#' @rdname runlist
#' @export
as.runlist <- function(x, ...) {
  do.call("runlist", as.list(x))
}


#' @export
"[.runlist" <- function(x, i, ...) {
  structure(NextMethod("["), class = class(x))
}


#' @rdname runlist
#' @export
c.hydromad <- function(..., recursive = FALSE) {
  args <- list(...)

  ## If some are runlists, then call c.runlist instead
  is.runlist <- sapply(args, inherits, what = "runlist")
  if (any(is.runlist)) {
    return(c.runlist(...))
  }

  ## If all hydromad objects, equivalent to just calling runlist
  if (!all(sapply(args, inherits, what = "hydromad"))) stop("Expected all elements to be hydromad objects")
  runlist(...)
}

#' @rdname runlist
#' @export
c.runlist <- function(..., recursive = FALSE) {
  args <- list(...)

  ## If any are hydromad objects, convert them to unit runlists
  is.hydromad <- sapply(args, inherits, what = "hydromad")
  for (i in which(is.hydromad)) args[[i]] <- do.call(runlist, args[i])

  ## Preserve classes, because primitive c strips them
  ## We keep only the classes common to all runlists
  ## Other user-defined attributes are still lost
  classes <- Reduce(intersect, lapply(args, class))

  ## Call primitive c
  args <- lapply(args, unclass)
  names(args) <- NULL
  rval <- do.call("c", args)

  ## Make a valid runlist
  rval <- as.runlist(rval)

  ## Restore classes
  class(rval) <- classes

  rval
}

#' @rdname runlist
#' @export
coef.runlist <-
  function(object, ..., items = NULL) {
    summary(object, ..., FUN = coef, items = items)
  }

#' @rdname runlist
#' @export
summary.runlist <-
  function(object, ..., FUN = summary, items = NULL) {
    stopifnot(is.list(object))
    if (length(object) == 0) {
      return(NULL)
    }
    ## extract elements from summary which are single numbers
    cc <- lapply(object, function(x, ...) {
      tmp <- FUN(x, ...)
      if (is.null(items)) {
        tmp <- tmp[unlist(lapply(tmp, function(z) {
          is.numeric(z) && !is.matrix(z) &&
            (length(z) == 1)
        }))]
      } else {
        tmp <- tmp[items]
      }
      unlist(tmp)
    }, ...)
    ## pad out missing entries with NAs
    ## find the set of all names
    allnms <- unique(unlist(lapply(cc, names)))
    ans <- matrix(NA_real_,
      nrow = length(object),
      ncol = length(allnms),
      dimnames = list(names(object), allnms)
    )
    for (i in 1:NROW(ans)) {
      ans[i, names(cc[[i]])] <- cc[[i]]
    }
    ans <- as.data.frame(ans)
    class(ans) <- c("summary.runlist", class(ans))
    ans
  }

# print.summary.runlist <-
#    function(x, digits = max(4, getOption("digits") - 3), ...)
# {
#    ## just simplify the printed output by rounding
#    print.data.frame(x, digits = digits, ...)
#    invisible(x)
# }

#' @rdname runlist
#' @export
print.runlist <-
  function(x, ...) {
    cat("\nList of model runs:\n")
    print.default(x, ...)
    invisible(x)
  }

#' @rdname runlist
#' @export
residuals.runlist <-
  function(object, ...) {
    ans <- lapply(object, residuals, ...)
    bad <- sapply(ans, length) == 0
    if (any(bad)) {
      stop(
        "residuals() returned nothing for items ",
        toString(names(ans)[bad])
      )
    }
    do.call("cbind", ans)
  }

#' @rdname runlist
#' @export
fitted.runlist <-
  function(object, ...) {
    ans <- lapply(object, fitted, ...)
    bad <- sapply(ans, length) == 0
    if (any(bad)) {
      stop(
        "fitted() returned nothing for items ",
        toString(names(ans)[bad])
      )
    }
    do.call("cbind", ans)
  }

#' @rdname runlist
#' @export
update.runlist <-
  function(object, ...) {
    switch(hydromad.getOption("parallel")[["update.runlist"]],
      "clusterApply" = {
        runs <- as.runlist(parLapply(cl, object, update, ...))
      },
      runs <- as.runlist(lapply(object, update, ...))
    ) ## switch parallel
    return(runs)
  }

#' @rdname runlist
#' @export
isValidModel.runlist <- function(object, ...) {
  return(TRUE)
}

#' @import utils
utils::globalVariables(c("cl"))
