## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

.defaultHydromadOptions <- function()
    list(
         sma = NULL,
         routing = NULL,
         scalar = scalar.ranges(),
         cwi = cwi.ranges(),
         cmd = cmd.ranges(),
         sacramento = sacramento.ranges(),
         bucket = bucket.ranges(),
         gr4j = gr4j.ranges(),
         gr4jrouting = gr4jrouting.ranges(),
         awbm = awbm.ranges(),
         snow = snow.ranges(),
         runoffratio = runoffratio.ranges(),
         dbm = dbm.ranges(),
         powuh = powuh.ranges(),
         order = c(n = 1, m = 0),
         delay = NA,
         max.delay = 10,
         warmup = 100,
         normalise = TRUE,
         loglik = function(Q, X, ...) -0.5 * sum((Q-X)^2, na.rm = TRUE),
         objective = ~ 0.7 * hmadstat("r.sq.sqrt")(Q, X) + 0.3 * .(hmadstat("r.sq.monthly", DATA = DATA))(Q, X),
         summary.stats = c("rel.bias", "r.squared", "r.sq.sqrt", "r.sq.log"),
         prefilter = list(0.9, 0.98, 0.5),
         sriv.iterations = 12,
         sriv.epsilon = 1e-3,
         riv.noise.order = NULL,
         inverse.fit.method = "sriv",
         inverse.iterations = 30,
         inverse.rel.tolerance = 1/5000,
         inverse.par.epsilon = 0.001,
         sce.control = list(),
         de.control = list(itermax = 1000/50),
         dream.control = list(),
         cmaes.control=list(),
         nsga2.control=list(),
         dds.control=list(
           logfile=NULL,
           projectfile=NULL,
           load_projectfile="no"
           ),
         fit.samples = 100,
         optim.method = "PORT",
         optim.control = list(reltol = 1e-6, maxit = 150,
                trace = if (interactive()) 1 else 0, REPORT = 10),
         nlminb.control = list(eval.max = 200, iter.max = 150, abs.tol = 0,
                trace = if (interactive()) 10 else 0, rel.tol = 1e-6),
         sim.epsilon = 1e-5,
         quiet = FALSE,
         trace = FALSE,
         catch.errors = TRUE,
         catch.errors.optim = TRUE,
         pure.R.code = FALSE,
	   parallel=list(update.runlist="none",objFunVal.runlist="none",evalPars="none",paretoObjectivesVaryWeights="none")
         )

## code below copied from lattice

hydromad.getOption <- function(name)
{
    .HydromadEnv$options[[name]]
}

hydromad.options <- function(...)
{
    ## this would have been really simple if only form allowed were
    ## lattice.options("foo", "bar") and
    ## lattice.options(foo=1, bar=2). But it could also be
    ## lattice.options(foo=1, "bar"), which makes some juggling necessary

    new <- list(...)
    if (is.null(names(new)) && length(new) == 1 && is.list(new[[1]])) new <- new[[1]]
    old <- .HydromadEnv$options

    ## if no args supplied, returns full options list
    if (length(new) == 0) return(old)

    nm <- names(new)
    if (is.null(nm)) return(old[unlist(new)]) ## typically getting options, not setting
    isNamed <- nm != "" ## typically all named when setting, but could have mix
    if (any(!isNamed)) nm[!isNamed] <- unlist(new[!isNamed])

    ## so now everything has non-"" names, but only the isNamed ones should be set
    ## everything should be returned, however

    retVal <- old[nm]
    names(retVal) <- nm
    nm <- nm[isNamed]

    ## this used to be

    ## modified <- updateList(retVal[nm], new[nm])
    ## .LatticeEnv$lattice.options[names(modified)] <- modified

    ## but then calling lattice.options(foo = NULL) had no effect
    ## because foo would be missing from modified.  So, we now do:

    updateList <- function (x, val) {
        if (is.null(x)) x <- list()
        utils::modifyList(x, val)
    }
    .HydromadEnv$options <- updateList(old, new[nm])

    ## return changed entries invisibly
    invisible(retVal)
}
