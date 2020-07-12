## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## a = time to drop to half
powuh.sim <-
  function(U, delay = 0,
           a = 1, b = 1, c = 1, init = 0, uhsteps = 100,
           na.action = na.pass,
           epsilon = hydromad.getOption("sim.epsilon")) {
    delay <- round(delay)
    ## note U is allowed to be multi-variate, i.e. multiple columns
    U <- na.action(U)
    ## apply 'U' delay in reverse to 'X' (i.e. lag by 'delay' steps)
    ## so can take delay as 0 for simulation purposes
    if (delay != 0) {
      U <- lag(U, -delay)
    }

    t <- seq.int(0, uhsteps - 1)
    uh <- (1 / (1 + (t / a)^(b / c)))^c
    ## normalise:
    uh <- uh / sum(uh)
    ## initialisation
    init <- matrix(init, nrow = uhsteps, ncol = NCOL(U))
    Upad <- rbind(init, as.matrix(U))
    Xpad <- filter(Upad, uh, sides = 1)
    # X <- U
    # coredata(X) <- coredata(filter(coredata(U), uh, sides = 1))
    X <- U
    X[] <- tail(Xpad, -uhsteps)

    ## align results to original input
    X <- shiftWindow(X, delay)

    ## zap simulated values smaller than epsilon
    X[X < epsilon] <- 0
    X
  }

powuh.ranges <- function() {
  list(
    a = c(0.01, 60),
    b = c(0.5, 3),
    c = c(0.5, 2)
  )
}

ssg.powuh <- function(theta) {
  1
}

normalise.powuh <- function(theta) {
  theta
}
