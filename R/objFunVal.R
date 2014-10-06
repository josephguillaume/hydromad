## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


objFunVal <- function(x, objective, ...)
    UseMethod("objFunVal")

objFunVal.default <-
    function(x, objective = hydromad.getOption("objective"),
             ..., nan.ok = FALSE)
{
    stopifnot(is.numeric(x) || is.data.frame(x))
    stopifnot(length(colnames(x)) > 0)
    ## these can be referred to in `objective`
    DATA <- x
    X <- x[,"X"]
    delayedAssign("Q", x[,"Q"])
    delayedAssign("U", x[,"U"])
    ## catch the .() function (used for cacheing, see hydromad.stats)
    ## normally it would not get through to here; evaluated by fitBy*() etc.
    ## But this may be needed if objFunVal() is called directly.
    assign(".", function(x) x)
    objFunVal1 <- function(obj, ...)
    {
        if (inherits(obj, "formula")) {
            val <- eval(obj[[2]])
        } else if (is.function(obj)) {
            assign(".", function(x) x, environment(obj))
            val <- obj(Q, X, ..., U = U, DATA = DATA)
        } else {
            stop("'objective' should be a function or formula, not a ",
                     toString(class(obj)))
        }
        if (is.nan(val)) {
            if (identical(nan.ok, "warn"))
                warning("objective function returned NaN")
            else if (!isTRUE(nan.ok))
                stop("objective function returned NaN")
        }
        if (!is.numeric(val))
            stop("objective should be numeric, not ", toString(class(val)))
        if (length(val) != 1)
            stop("objective value should be length 1, not ", length(val))
        as.numeric(val)
    }
    if (is.list(objective))
        lapply(objective, objFunVal1, ...)
    else
        objFunVal1(objective, ...)
}

## TODO: could this just merge the data and call the default method? slow?
## 
objFunVal.tf <-
objFunVal.hydromad <-
    function(x, objective = hydromad.getOption("objective"),
             ..., all = FALSE, nan.ok = FALSE)
{
    model <- x
    ## these can be referred to in `objective`
    X <- fitted(x, all = all)
    if (length(X) == 0)
        stop("fitted() returned nothing")
    delayedAssign("Q", observed(x, all = all))
    delayedAssign("U", fitted(x, all = all, U = TRUE))
    delayedAssign("DATA", observed(x, all = all, select = TRUE))
    ## catch the .() function (used for cacheing, see hydromad.stats)
    ## normally it would not get through to here; evaluated by fitBy*() etc
    ## But this may be needed if objFunVal() is called directly.
    assign(".", function(x) x)
    isValidModel <- isValidModel(x)
    objFunVal1 <- function(obj, ...)
    {
        if (!isValidModel)
            return(NA_real_)
        if (inherits(obj, "formula")) {
            val <- eval(obj[[2]])
        } else if (is.function(obj)) {
            assign(".", function(x) x, environment(obj))
            val <- obj(Q, X, ..., U = U, DATA = DATA, model = model)
        } else {
            stop("'objective' should be a function or formula, not a ",
                 toString(class(obj)))
        }
        if (is.nan(val)) {
            if (identical(nan.ok, "warn"))
                warning("objective function returned NaN")
            else if (!isTRUE(nan.ok))
                stop("objective function returned NaN")
        }
        if (!is.numeric(val))
            stop("objective should be numeric, not ", toString(class(val)))
        if (length(val) != 1)
            stop("objective value should be length 1, not ", length(val))
        as.numeric(val)
    }
    if (is.list(objective))
        lapply(objective, objFunVal1, ...)
    else
        objFunVal1(objective, ...)
}

objFunVal.runlist <- function (x, objective = list(hydromad.getOption("objective"),mean), ...) 
{
  if(is.list(objective) && length(objective) == 2){
    switch(hydromad.getOption("parallel")[["objFunVal.runlist"]],
           clusterApply = {
             vals <- parSapply(cl,x, objFunVal, objective[[1]], ...)
           },
           vals <- sapply(x, objFunVal, objective[[1]], ...)
           )
    agg <- objective[[2]](vals)
    stopifnot(length(agg) == 1)
    return(agg)
  } else if(is.function(objective)){
    return(objective(x,...))
  }
  stop("Objective is not a list of length 2 or a function")
}
