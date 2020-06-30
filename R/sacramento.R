## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

## Sacramento Soil Moisture Accounting model.
## Developed by the US National Weather Service.
sacramento.sim <-
    function(DATA,
             uztwm,
             uzfwm,
             uzk,
             pctim,
             adimp,
             zperc,
             rexp,
             lztwm,
             lzfsm,
             lzfpm,
             lzsk,
             lzpk,
             pfree,
             etmult = 1,
             dt = 1,
             uztwc_0 = 0.5,
             uzfwc_0 = 0.5,
             lztwc_0 = 0.5,
             lzfsc_0 = 0.5,
             lzfpc_0 = 0.5,
             adimc_0 = 0.5,
             min_ninc = 20,
             return_state = FALSE)
    {
        stopifnot(c("P", "E") %in% colnames(DATA))
        ## check values
        stopifnot(uztwm >= 0)
        stopifnot(uzfwm >= 0)
        stopifnot(uzk >= 0)
        stopifnot(0 <= pctim && pctim <= 1)
        stopifnot(adimp >= 0)
        stopifnot(zperc >= 0)
        stopifnot(lztwm >= 0)
        stopifnot(lzfsm >= 0)
        stopifnot(lzfpm >= 0)
        stopifnot(lzsk >= 0)
        stopifnot(lzpk >= 0)
        stopifnot(pfree >= 0)
        stopifnot(etmult >= 0)
        stopifnot(dt >= 0)
        
        xpar <-
            c(uztwm,
              uzfwm,
              uzk,
              pctim,
              adimp,
              zperc,
              rexp,
              lztwm,
              lzfsm,
              lzfpm,
              lzsk,
              lzpk,
              pfree)
        
        P <- DATA[, "P"]
        E <- DATA[, "E"]
        ## skip over missing values
        bad <- is.na(P) | is.na(E)
        P[bad] <- 0
        E[bad] <- 0
        ## NOTE: there is no pure.R.code implementation
        COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
        if (COMPILED) {
            if (return_state) {
                states <- .C(
                    sma_sac_state,
                    as.double(P),
                    as.double(E),
                    as.integer(NROW(DATA)),
                    as.double(xpar),
                    as.double(etmult),
                    as.double(dt),
                    U = double(NROW(DATA)),
                    uztwc = double(NROW(DATA)),
                    uzfwc = double(NROW(DATA)),
                    lztwc = double(NROW(DATA)),
                    lzfsc = double(NROW(DATA)),
                    lzfpc = double(NROW(DATA)),
                    adimc = double(NROW(DATA)),
                    sett = double(NROW(DATA)),
                    se1 = double(NROW(DATA)),
                    se3 = double(NROW(DATA)),
                    se4 = double(NROW(DATA)),
                    se5 = double(NROW(DATA)),
                    roimp = double(NROW(DATA)),
                    sdro = double(NROW(DATA)),
                    ssur = double(NROW(DATA)),
                    sif = double(NROW(DATA)),
                    bfp = double(NROW(DATA)),
                    bfs = double(NROW(DATA)),
                    bfcc = double(NROW(DATA)),
                    as.double(uztwc_0),
                    as.double(uzfwc_0),
                    as.double(lztwc_0),
                    as.double(lzfsc_0),
                    as.double(lzfpc_0),
                    as.double(adimc_0),
                    as.integer(min_ninc),
                    NAOK = FALSE,
                    DUP = FALSE,
                    PACKAGE = "hydromad"
                )
                for (i in 7:25)
                    attributes(states[[i]]) <- attributes(P)
                ans <- do.call(cbind, states[7:25])
                return(ans)
            } else {
                U <- .C(
                    sma_sac,
                    as.double(P),
                    as.double(E),
                    as.integer(NROW(DATA)),
                    as.double(xpar),
                    as.double(etmult),
                    as.double(dt),
                    U = double(NROW(DATA)),
                    as.double(uztwc_0),
                    as.double(uzfwc_0),
                    as.double(lztwc_0),
                    as.double(lzfsc_0),
                    as.double(lzfpc_0),
                    as.double(adimc_0),
                    as.integer(min_ninc),
                    NAOK = FALSE,
                    DUP = FALSE,
                    PACKAGE = "hydromad"
                )$U
                ## make it a time series object again
                attributes(U) <- attributes(P)
                ## re-insert missing values
                U[bad] <- NA
                return(U)
            }
        } else {
            U <- rep(0, length(P))
            sma = list(
                uztwm = uztwm,
                uzfwm = uzfwm,
                uzk = uzk,
                pctim = pctim,
                adimp = adimp,
                riva  = 0.0,
                zperc = zperc,
                rexp  = rexp,
                lztwm = lztwm,
                lzfsm = lzfsm,
                lzfpm = lzfpm,
                lzsk  = lzsk,
                lzpk  = lzpk,
                pfree = pfree,
                rserv = 0.3,
                side  = 0.0,
                pxmlt = 1.0,
                uztwc = uztwc_0 * uztwm,
                uzfwc = uzfwc_0 * uzfwm,
                lztwc = lztwc_0 * lztwm,
                lzfsc = lzfsc_0 * lzfsm,
                lzfpc = lzfpc_0 * lzfpm,
                adimc = adimc_0 * (uztwm + lztwm),
                min_ninc = min_ninc
            )
            
            # SET SOME INITIAL VALUES TO ZERO
            fsum1 = list(
                srot = 0,
                simpvt = 0,
                srodt = 0,
                srost = 0,
                sintft = 0,
                sgwfp = 0,
                sgwfs = 0,
                srecht = 0,
                sett = 0,
                se1 = 0,
                se3 = 0,
                se4 = 0,
                se5 = 0
            )
            sma$epdist = etmult
            sma$dt = dt
            for (t in 1:NROW(DATA)) {
                sma$ep =  E[t] * sma$epdist
                sma$pxv = P[t] * sma$pxmlt
                fland1(sma,fsum1)
                U[t] = sma$tlci
                t<-t+1
            }
            # TODO: Construct sma and fsum1 as in test_fland1.R
            # TODO: for loop as in sacramento.c
            #    TODO: inside for loop, call fland1
            #    TODO:  set the value of U for that timestep
            
            
            ## make it a time series object again
            attributes(U) <- attributes(P)
            ## re-insert missing values
            U[bad] <- NA
            return(U)
        }
    }

sacramento.ranges <- function()
    list(
        uztwm = c(1, 150),
        uzfwm = c(1, 150),
        uzk = c(0.1, 0.5),
        pctim = c(0.000001, 0.1),
        adimp = c(0, 0.4),
        zperc = c(1, 250),
        rexp = c(0, 5),
        lztwm = c(1, 500),
        lzfsm = c(1, 1000),
        lzfpm = c(1, 1000),
        lzsk = c(0.01, 0.25),
        lzpk = c(0.0001, 0.25),
        pfree = c(0, 0.6)
    )