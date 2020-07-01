## SWIMP: Simple Wetland Inundation Model using Poweroids
## by Felix Andrews, Barry Croke and Baihua Fu 2009
## (inspired by MFAT Floodplain Hydrology model)

swimp <-
  function(flow.ML,
           thresh = 0,
           alpha,
           beta,
           E.mm = 0,
           P.mm = 0,
           Ksat.mm.day = 0,
           e = 0.2,
           g = 140,
           Hmax = 2000,
           Amax = 10000,
           porosity = 0.2,
           M_0 = Hmax * porosity,
           V_0 = 0,
           drainage = 0,
           drainLevel = 0) {
    stopifnot(NROW(flow.ML) > 1)
    stopifnot(NCOL(flow.ML) == 1)
    force(list(alpha, beta))
    overflow <- pmax(coredata(flow.ML) - thresh, 0)
    E.mm <- rep(coredata(E.mm), length = length(overflow))
    P.mm <- rep(coredata(P.mm), length = length(overflow))
    ## volume of surface water, ML
    V <- double(length(overflow))
    ## area of surface water, km^2
    A <- double(length(overflow))
    ## surface water level, mm
    H <- double(length(overflow))
    ## ground water level as moisture deficit (CMD), mm
    M <- double(length(overflow))
    ## infiltration into wetland
    Iw <- double(length(overflow))
    ## initial conditions
    M[1] <- M_0
    V[1] <- V_0
    ## work out volume below drainage level
    undrainedV <- poweroid(H = drainLevel, alpha = alpha, beta = beta)$V

    ## TODO: start from t=1 (add dummy row to time series?)
    if (isTRUE(hydromad.getOption("pure.code"))) {
      ## slower version in R for cross-checking
      storedV <- poweroid(H = drainLevel, alpha = alpha, beta = beta)$V
      for (t in 2:NROW(overflow)) {
        ## work out water level in ground
        ## relative to base of wetland, H==0
        Hg <- Hmax - (1 / porosity) * M[t - 1]
        ## infiltration rate (drainage) from wetland
        ## depends on pressure; reference = 100mm
        Iw[t] <- Ksat.mm.day * (H[t - 1] - Hg) / 100
        ## if discharging into wetland, adjust mass for porosity
        if (H[t - 1] < Hg) Iw[t] <- Iw[t] * porosity
        ## evapo-transpiration (equation from CMD model)
        Mf <- M[t - 1] - P.mm[t]
        ETg <- e * E.mm[t] * min(1, exp(2 * (1 - Mf / g)))
        ETg <- max(ETg, 0)
        ## mass balance of water level in ground (CMD)
        M[t] <- M[t - 1] + ETg - P.mm[t] - Iw[t] * (A[t - 1] / Amax)
        M[t] <- max(M[t], 0)
        ## mass balance of volume in wetland
        V[t] <- V[t - 1] + overflow[t] + (P.mm[t] - E.mm[t]) * A[t - 1] - Iw[t] * A[t - 1]
        V[t] <- max(V[t], 0)
        ## drainage of volume above drainage level
        V[t] <- V[t] - max(V[t] - undrainedV, 0) * drainage
        ## convert volume to water level and area
        ## NOTE: these formulae are designed to avoid numerical explosions!
        ## e.g. ^(1/beta) when beta = 0.01
        H[t] <- ((2 / beta + 1) * V[t] / pi)^(beta / (2 + beta)) * alpha^(2 / (2 + beta))
        A[t] <- ((2 / beta + 1) * V[t] / (alpha * pi^(-beta / 2)))^(2 / (2 + beta))
      }
    } else {
      ans <- .C(swimp_core,
        as.double(overflow),
        as.integer(length(overflow)),
        as.double(alpha),
        as.double(beta),
        as.double(E.mm),
        as.double(P.mm),
        as.double(Ksat.mm.day),
        as.double(e),
        as.double(g),
        as.double(Hmax),
        as.double(Amax),
        as.double(porosity),
        as.double(drainage),
        as.double(undrainedV),
        V = V,
        A = A,
        H = H,
        M = M,
        Iw = Iw,
        DUP = FALSE, PACKAGE = "hydromad"
      )
      V <- ans$V
      A <- ans$A
      H <- ans$H
      M <- ans$M
      Iw <- ans$Iw
    }
    mean.depth <- ifelse(A > 0, V / A, 0)
    ans <- zoo(
      cbind(
        volume = V, area = A,
        level = H, mean.depth = mean.depth
      ),
      time(flow.ML)
    )
    if (Ksat.mm.day > 0) {
      ans <- cbind(ans, Iw = Iw, CMD = M)
    }
    ans
  }
