## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## local environment to store user options
.HydromadEnv <- new.env()
.HydromadEnv$options <- list()
.HydromadEnv$stats <- list()

.onLoad <- function(libname, pkgname) {
  hydromad.options(.defaultHydromadOptions())
  hydromad.stats(.defaultHydromadStats())
}
