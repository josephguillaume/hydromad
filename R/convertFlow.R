##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

convertFlow <-
  function(x, from = "mm", to = "mm", area.km2 = -1,
           timestep.default = "days") {
    if (from == "cumecs") from <- "m^3/sec"
    if (to == "cumecs") to <- "m^3/sec"

    ## extract timestep from each of 'from' and 'to'
    from.step <- to.step <- timestep.default
    if (any(grep("/", from))) {
      from.step <- sub("^.*/ *", "", from)
      from <- sub(" */.*$", "", from)
    }
    if (any(grep("/", to))) {
      to.step <- sub("^.*/ *", "", to)
      to <- sub(" */.*$", "", to)
    }
    ## extract multiplier from each timestep
    from.mult <- gsub("[^0-9\\.]", "", from.step)
    from.step <- gsub("[0-9\\. ]", "", from.step)
    to.mult <- gsub("[^0-9\\.]", "", to.step)
    to.step <- gsub("[0-9\\. ]", "", to.step)
    ## number of seconds for each possible time step
    timefactors <- alist(
      millisecond = , ms = 0.001,
      seconds = , second = , sec = , s = 1,
      minutes = , minute = , min = 60,
      hours = , hour = , hr = , h = 60 * 60,
      days = , day = , d = 24 * 60 * 60,
      weeks = , week = 7 * 24 * 60 * 60,
      months = , month = , mon = 30.4375 * 24 * 60 * 60,
      annum = , anna = , a = , years = , year = , yr = , y = 365.25 * 24 * 60 * 60
    )
    from.secs <- do.call(switch, c(from.step, timefactors))
    to.secs <- do.call(switch, c(to.step, timefactors))
    if (is.null(from.secs) || is.null(to.secs)) {
      stop("unrecognised time unit")
    }
    if (nchar(from.mult) > 0) from.secs <- from.secs * as.numeric(from.mult)
    if (nchar(to.mult) > 0) to.secs <- to.secs * as.numeric(to.mult)
    ## handle volumes
    depthUnits <- c("mm", "cm", "metres", "km", "inches", "feet", "ft")
    volUnits <- c(
      "mL", "cL", "dL", "L", "daL", "hL", "kL", "ML", "GL", "TL",
      "cm3", "dm3", "m3", "km3", "ft3",
      "cm^3", "dm^3", "m^3", "km^3", "ft^3"
    )
    allUnits <- c(depthUnits, volUnits)
    from <- match.arg(from, allUnits)
    to <- match.arg(to, allUnits)
    if ((from %in% depthUnits) != (to %in% depthUnits)) {
      if (missing(area.km2)) stop("need to give 'area.km2'")
    }

    ## factors to convert to mm (*) or from mm (/) per timestep
    Litres <- (1 / area.km2) / 1e6
    vfactors <- alist(
      mm = 1,
      cm = 10,
      metres = , metre = , m = 1000,
      km = 1000 * 1000,
      inches = , inch = , `in` = 25.4,
      feet = , ft = 304.8,
      mL = , cm3 = , `cm^3` = 1e-3 * Litres,
      cL = 0.01 * Litres,
      dL = 0.1 * Litres,
      L = , dm3 = , `dm^3` = Litres,
      daL = 10 * Litres,
      hL = 100 * Litres,
      kL = , m3 = , `m^3` = 1000 * Litres,
      ML = 1e6 * Litres,
      GL = 1e9 * Litres,
      TL = , km3 = , `km^3` = 1e12 * Litres,
      ft3 = , `ft^3` = 1 / 0.0353146667 * Litres,
      stop("unrecognised volume unit")
    )
    ## first convert to mm
    x <- x * do.call(switch, c(from, vfactors))
    ## now convert to required unit 'to'
    x <- x / do.call(switch, c(to, vfactors))
    ## now convert timesteps
    x <- x * (to.secs / from.secs)
    x
  }
