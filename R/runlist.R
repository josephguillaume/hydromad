## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

runlist <- function(...)
{
    object <- list(...)
    if (is.null(names(object)))
        names(object) <- rep("", length(object))
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
    if (length(object) > 0)
        class(object) <- paste(class(object[[1]]), "runlist", sep = ".")
    class(object) <- unique(c(class(object), "runlist", "list"))
    object
}

as.runlist <- function(x, ...)
    do.call("runlist", as.list(x))

"[.runlist" <- function (x, i, ...)
    structure(NextMethod("["), class = class(x))

c.runlist <- function(..., recursive = FALSE)
    as.runlist(NextMethod("c"))

coef.runlist <-
    function(object, ..., items = NULL)
{
    summary(object, ..., FUN = coef, items = items)
}


summary.runlist <-
    function(object, ..., FUN = summary, items = NULL)
{
    stopifnot(is.list(object))
    if (length(object) == 0)
        return(NULL)
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
                  dimnames = list(names(object), allnms))
    for (i in 1:NROW(ans))
        ans[i, names(cc[[i]])] <- cc[[i]]
    ans <- as.data.frame(ans)
    class(ans) <- c("summary.runlist", class(ans))
    ans
}

#print.summary.runlist <-
#    function(x, digits = max(4, getOption("digits") - 3), ...)
#{
#    ## just simplify the printed output by rounding
#    print.data.frame(x, digits = digits, ...)
#    invisible(x)
#}

print.runlist <-
    function(x, ...)
{
    cat("\nList of model runs:\n")
    print.default(x, ...)
    invisible(x)
}

residuals.runlist <-
    function(object, ...)
{
    ans <- lapply(object, residuals, ...)
    bad <- sapply(ans, length) == 0
    if (any(bad))
        stop("residuals() returned nothing for items ",
             toString(names(ans)[bad]))
    do.call("cbind", ans)
}

fitted.runlist <-
    function(object, ...)
{
    ans <- lapply(object, fitted, ...)
    bad <- sapply(ans, length) == 0
    if (any(bad))
        stop("fitted() returned nothing for items ",
             toString(names(ans)[bad]))
    do.call("cbind", ans)
}

update.runlist <-
    function(object, ...)
{
  switch(hydromad.getOption("parallel")[["update.runlist"]],
         "clusterApply"={
           runs <- as.runlist(parLapply(cl,object,update,...))
         },
         runs <- as.runlist(lapply(object, update, ...))
         ) ## switch parallel
  return(runs)
}

isValidModel.runlist <- function(object,...)
  return(TRUE)
