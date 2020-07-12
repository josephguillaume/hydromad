fland1 <- function(sma, fsum1) {
  EPS <- 1e-5
  zp <- sma$zperc  #COMPUTE EVAPOTRANSPIRATION LOSS FOR THE TIME INTERVAL
  #EDMND IS THE ET-DEMAND FOR THE TIME INTERVAL

  edmnd <- sma$ep
  e1 <- edmnd * sma$uztwc / sma$uztwm
  red <- edmnd - e1 #red is residual evap demand
  sma$uztwc <- sma$uztwc - e1
  if (abs(sma$uztwc) < EPS) {
    sma$uztwc <- 0
  }
  e2 <- 0
  #e1 cannot exceed sma$uztwc
  if (sma$uztwc < 0) {
    e1 <- e1 + sma$uztwc
    sma$uztwc <- 0
    red <- edmnd - e1
    if (sma$uzfwc < red) {
      e2 <- sma$uzfwc  #e2 is evap from sma$uzfwc
      sma$uzfwc <- 0
      red <- red - e2
    }
    else {
      e2 <- red
      sma$uzfwc <- sma$uzfwc - e2
      red <- 0
      #UPPER ZONE FREE WATER RATIO EXCEEDS UPPER ZONE TENSION WATER RATIO, THUS TRANSFER FREE WATER TO TENSION
      if ((sma$uztwc / sma$uztwm) < (sma$uzfwc / sma$uzfwm)) {
        uzrat <- (sma$uztwc + sma$uzfwc) / (sma$uztwm + sma$uzfwm)
        sma$uztwc <- sma$uztwm * uzrat
        sma$uzfwc <- sma$uzfwm * uzrat
      }
    }
  }

  # L225
  e3 <- red * sma$lztwc / (sma$uztwm + sma$lztwm)
  sma$lztwc <- sma$lztwc - e3
  if (abs(sma$lztwc) < EPS) {
    sma$lztwc <- 0
  }
  if (sma$lztwc < 0) {
    e3 <- e3 + sma$lztwc #e3 cannot exceed sma$lztwc
    sma$lztwc <- 0
  }
  ratlzt <- sma$lztwc / sma$lztwm
  saved <- sma$rserv * (sma$lzfpm + sma$lzfsm)
  ratlz <- (sma$lztwc + sma$lzfpc + sma$lzfsc - saved) / (sma$lztwm + sma$lzfpm +
    sma$lzfsm - saved)

  if (ratlzt < ratlz) {
    #RESUPPLY LOWER ZONE TENSION WATER FROM LOWER ZONE FREE WATER IF MORE WATER AVAILABLE THERE
    del <- (ratlz - ratlzt) * sma$lztwm
    #TRANSFER FROM sma$lzfsc to sma$lztwc
    sma$lztwc <- sma$lztwc + del
    sma$lzfsc <- sma$lzfsc - del
    if (sma$lzfsc < 0) {  #IF TRANSFER EXCEEDS LZFSC THEN REMAINDER COMES FROM LZFPC
      sma$lzfpc <- sma$lzpfpc + sma$lzfsc
      sma$lzfsc <- 0
    }
  }
  #COMPUTE ET FROM ADIMP AREA.-E5
  e5 <- e1 + (red + e2) * (sma$adimc - e1 - sma$uztwc) / (sma$uztwm + sma$lztwm)
  
  #ADJUST ADIMC,ADDITIONAL IMPERVIOUS AREA STORAGE, FOR EVAPORATION
  sma$adimc <- sma$adimc - e5
  if (abs(sma$adimc) < EPS) {
    sma$adimc <- 0
  }
  if (sma$adimc < 0) {  #E5 CAN NOT EXCEED ADIMC
    e5 <- e5 + sma$adimc
    sma$adimc <- 0
  }
  e5 <- e5 * sma$adimp
  #E5 IS ET FROM THE AREA ADIMP.
  #COMPUTE PERCOLATION AND RUNOFF AMOUNTS

  twx <- sma$pxv + sma$uztwc - sma$uztwm
  if (twx < 0) {  #ALL MOISTURE HELD IN UZTW--NO EXCESS
    sma$uztwc <- sma$uztwc + sma$pxv
    twx <- 0
  }
  else { #MOISTURE AVAILABLE IN EXCESS OF UZTW STORAGE
    sma$uztwc <- sma$uztwm
  }
  sma$adimc <- sma$adimc + (sma$pxv - twx)
  #COMPUTE IMPERVIOUS AREA RUNOFF
  roimp <- sma$pxv * sma$pctim
  #ROIMP IS RUNOFF FROM THE MINIMUM IMPERVIOUS AREA
  fsum1$simpvt <- fsum1$simpvt - roimp
  
  #INITIALIZE TIME INTERVAL SUMS
  sbf <- ssur <- sif <- sperc <- sdro <- spbf <- 0
  
  #DETERMINE COMPUTATIONAL TIME INCREMENTS FOR THE BASIC TIME INTERVAL
  # NINC = NUMBER OF TIME INCREMENTS THAT THE TIME INTERVAL
  # IS DIVIDED INTO FOR FURTHER SOIL-MOISTURE ACCOUNTING.
  # NO ONE INCREMENT WILL EXCEED 5.0 MILLIMETERS OF UZFWC+PAV
  # DINC = LENGTH OF EACH INCREMENT IN DAYS.
  # PINC = AMOUNT OF AVAILABLE MOISTURE FOR EACH INCREMENT.
  # ninc = (int) (1.0 + 0.2 * (sma->uzfwc + twx));

  ninc <- as.integer(1.0 + 0.2 * (sma$uzfwc + twx))
  if (ninc < sma$min_ninc) { #minimum number of loops
    ninc <- sma$min_ninc
  }
  dinc <- 1.0 / ninc * sma$dt
  pinc <- twx / ninc
  #COMPUTE FREE WATER DEPLETION FRACTIONS FOR THE TIME INCREMENT 
  #BEING USED-BASIC DEPLETIONS ARE FOR ONE DAY

  if (sma$uzk > 1) {
    print(sma$uzk)
  }
  if (sma$lzpk > 1) {
    print(sma$lzpk)
  }
  if (sma$lzsk > 1) {
    print(sma$lzsk)
  }
  duz <- 1.0 - (1.0 - sma$uzk)^dinc
  dlzp <- 1.0 - (1.0 - sma$lzpk)^dinc
  dlzs <- 1.0 - (1.0 - sma$lzsk)^dinc
  parea <- 1.0 - sma$adimp - sma$pctim

  #START INCREMENTAL FOR LOOP FOR THE TIME INTERVAL
  for (i in 1:ninc) {
    adsur <- 0
    
    #COMPUTE DIRECT RUNOFF (FROM ADIMP AREA).
    #ADDRO IS THE AMOUNT OF DIRECT RUNOFF FROM THE AREA ADIMP-SDRO IS THE SIX HOUR SUMMATION
    ratio <- (sma$adimc - sma$uztwc) / sma$lztwm
    addro <- pinc * (ratio * ratio)
    sdro <- sdro + (addro * sma$adimp)
    
    #COMPUTE BASEFLOW AND KEEP TRACK OF TIME INTERVAL SUM
    bf <- sma$lzfpc * dlzp
    sma$lzfpc <- sma$lzfpc - bf
    if (sma$lzfpc <= 1e-4) {
      bf <- bf + sma$lzfpc
      sma$lzfpc <- 0
    }
    sbf <- sbf + bf
    spbf <- spbf + bf
    bf <- sma$lzfsc * dlzs
    sma$lzfsc <- sma$lzfsc - bf
    if (sma$lzfsc <= 1e-4) {
      bf <- bf + sma$lzfsc
      sma$lzfsc <- 0
    }
    sbf <- sbf + bf
    
    #COMPUTE PERCOLATION-IF NO WATER AVAILABLE THEN SKIP
    if ((pinc + sma$uzfwc) <= 1e-2) {
      sma$uzfwc <- sma$uzfwc + pinc
      sma$adimc <- sma$adimc + (pinc - addro - adsur) # L249
      if (sma$adimc > (sma$uztwm + sma$lztwm)) {
        addro <- addro + (sma$adimc - (sma$uztwm + sma$lztwm))
        sma$adimc <- sma$uztwm + sma$lztwm
      }
      next
    }
    percm <- sma$lzfpm * dlzp + sma$lzfsm * dlzs

    if (zp < 0.0) {
      zp <- 0.0
    }
    perc <- percm * sma$uzfwc / sma$uzfwm
    
    #DEFR IS THE LOWER ZONE MOISTURE DEFICIENCY RATIO
    defr <- 1.0 - (sma$lztwc + sma$lzfpc + sma$lzfsc) / (sma$lztwm + sma$lzfpm + sma$lzfsm)

    if (defr < 0) {
      print(defr)
      c <- c(sma$lztwc, sma$lzfpc, sma$lzfsc)
      print(c)
      m <- c(sma$lztwm, sma$lzfpm, sma$lzfsm)
      print(m)
      stop("1") # doubt
    }

    perc <- perc * (1.0 + zp * (defr^sma$rexp))
    
    #NOTE...PERCOLATION OCCURS FROM UZFWC BEFORE PAV IS ADDED
    if (perc >= sma$uzfwc) {  #PERCOLATION RATE EXCEEDS UZFWC
      perc <- sma$uzfwc
    }
    sma$uzfwc <- sma$uzfwc - perc  #PERCOLATION RATE IS LESS THAN UZFWC

    check <- sma$lztwc + sma$lzfpc + sma$lzfsc + perc - sma$lztwm - sma$lzfpm -
      sma$lzfsm
    if (check > 0) {  #CHECK TO SEE IF PERCOLATION EXCEEDS LOWER ZONE DEFICIENCY
      perc <- perc - check
      sma$uzfwc <- sma$uzfwc + check
    }
    
    #SPERC IS THE TIME INTERVAL SUMMATION OF PERC
    sperc <- sperc + perc
    
    #COMPUTE INTERFLOW AND KEEP TRACK OF TIME INTERVAL SUM.
    # NOTE...PINC HAS NOT YET BEEN ADDED
    del <- sma$uzfwc * duz
    sif <- sif + del
    sma$uzfwc <- sma$uzfwc - del
    
    #DESCRIBE PERCOLATED WATER INTO THE LOWER ZONES
    #TENSION WATER MUST BE FILLED FIRST EXCEPT FOR THE
    #PFREE AREA.  PERCT IS PERCOLATION TO TENSION WATER
    #AND PERCF IS PERCOLATION GOING TO FREE WATER

    perct <- perc * (1.0 - sma$pfree)
    if ((perct + sma$lztwc) <= sma$lztwm) {
      sma$lztwc <- sma$lztwc + perct
      percf <- 0
    }
    else {
      percf <- perct + sma$lztwc - sma$lztwm
      sma$lztwc <- sma$lztwm
    }
    
    #DISTRIBUTE PERCOLATION IN EXCESS OF TENSION REQUIREMENTS AMONG THE FREE WATER STORAGES

    percf <- percf + (perc * sma$pfree)
    if (percf != 0) {
      #HPL IS THE RELATIVE SIZE OF THE PRIMARY STORAGE AS COMPARED WITH TOTAL LOWER ZONE FREE WATER STORAGE.
      hpl <- sma$lzfpm / (sma$lzfpm + sma$lzfsm)
      
      #RATLP AND RATLS ARE CONTENT TO CAPACITY RATIOS, OR
      #IN OTHER WORDS, THE RELATIVE FULLNESS OF EACH STORAGE
      ratlp <- sma$lzfpc / sma$lzfpm
      ratls <- sma$lzfsc / sma$lzfsm

      #FRACP IS THE FRACTION GOING TO PRIMARY
      fracp <- hpl * 2.0 * (1.0 - ratlp) / (1.0 - ratlp + 1.0 - ratls)
      if (fracp > 1.0) {
        fracp <- 1.0
      }
      
      #PERCP AND PERCS ARE THE AMOUNT OF THE EXCESS PERCOLATION GOING TO 
      #PRIMARY AND SUPPLEMENTAL STORAGES, RESPECTIVELY.
      percp <- percf * fracp
      percs <- percf - percp
      sma$lzfsc <- sma$lzfsc + percs
      if (sma$lzfsc > sma$lzfsm) {
        percs <- percs + (-sma$lzfsc + sma$lzfsm)
        sma$lzfsc <- sma$lzfsm
      }
      sma$lzfpc <- sma$lzfpc + (percf - percs)

      #CHECK TO MAKE SURE LZFPC DOES NOT EXCEED LZFPM
      if (sma$lzfpc > sma$lzfpm) {
        excess <- sma$lzfpc - sma$lzfpm
        sma$lztwc <- sma$lztwc + excess
        sma$lzfpc <- sma$lzfpm
      }
    }
    
    #DISTRIBUTE PINC BETWEEN UZFWC AND SURFACE RUNOFF
    if (pinc != 0) {
      if ((pinc + sma$uzfwc) <= sma$uzfwm) {  #CHECK IF PINC EXCEEDS UZFWM
        sma$uzfwc <- sma$uzfwc + pinc  #NO SUFACE RUNOFF
      }
      else {
        sur <- pinc + sma$uzfwc - sma$uzfwm
        sma$uzfwc <- sma$uzfwm
        
        #ADSUR IS THE AMOUNT OF SURFACE RUNOFF WHICH COMES FROM 
        #THAT PORTION OF ADIMP WHICH IS NOT CURRENTLY GENERATING DIRECT RUNOFF.
        #ADDRO/PINC IS THE FRACTION OF ADIMP CURRENTLY GENERATING DIRECT RUNOFF.
        ssur <- ssur + (sur * parea)
        adsur <- sur * (1.0 - addro / pinc)
        ssur <- ssur + (adsur * sma$adimp)
      }
    }
    sma$adimc <- sma$adimc + (pinc - addro - adsur) # L249
    if (sma$adimc > (sma$uztwm + sma$lztwm)) {
      addro <- addro + (sma$adimc - (sma$uztwm + sma$lztwm))
      sma$adimc <- sma$uztwm + sma$lztwm
    }
  }  #END OF INCREMENTAL FOR LOOP
  
  #COMPUTE SUMS AND ADJUST RUNOFF AMOUNTS BY THE AREA OVER WHICH THEY ARE GENERATED.
  
  #EUSED IS THE ET FROM PAREA WHICH IS 1.0-ADIMP-PCTIM
  eused <- e1 + e2 + e3
  sif <- sif * parea
  
  #SEPARATE CHANNEL COMPONENT OF BASEFLOW FROM THE NON-CHANNEL COMPONENT
  tbf <- sbf * parea   #TBF IS TOTAL BASEFLOW
  
  #BFCC IS BASEFLOW, CHANNEL COMPONENT
  bfcc <- tbf * (1.0 / (1.0 + sma$side))
  bfp <- (spbf * parea) / (1.0 + sma$side)
  bfs <- bfcc - bfp
  if (bfs < 0) {
    bfs <- 0
  }
  bfncc <- tbf - bfcc  #BFNCC IS BASEFLOW, NON-CHANNEL COMPONENT
  
  #ADD TO MONTHLY SUMS
  fsum1$sintft <- fsum1$sintft + sif
  fsum1$sgwfp <- fsum1$sgwfp + bfp
  fsum1$sgwfs <- fsum1$sgwfs + bfs
  fsum1$srecht <- fsum1$srecht + bfncc
  fsum1$srost <- fsum1$srost + ssur
  fsum1$srodt <- fsum1$srodt + sdro

  #STORE EACH OF THE FIVE FLOWS IN VARIABLE TLCI_FLOWS
  sma$tlci_flows[1] <- roimp
  sma$tlci_flows[2] <- sdro
  sma$tlci_flows[3] <- ssur
  sma$tlci_flows[4] <- sif
  sma$tlci_flows[5] <- bfp
  sma$tlci_flows[6] <- bfs
  sma$tlci_flows[7] <- bfcc

  #COMPUTE TOTAL CHANNEL INFLOW FOR THE TIME INTERVAL
  sma$tlci <- roimp + sdro + ssur + sif + bfcc

  #COMPUTE E4-ET FROM RIPARIAN VEGETATION
  e4 <- (edmnd - eused) * sma$riva

  #SUBTRACT E4 FROM CHANNEL INFLOW
  sma$tlci <- sma$tlci - e4

  if (sma$tlci < 0) {
    e4 <- e4 + sma$tlci
    sma$tlci <- 0
  }
  fsum1$srot <- fsum1$srot + sma$tlci

  #COMPUTE TOTAL EVAPOTRANSPIRATION-TET
  eused <- eused * parea
  tet <- eused + e5 + e4
  fsum1$sett <- fsum1$sett + tet
  fsum1$se1 <- fsum1$se1 + (e1 * parea)
  fsum1$se3 <- fsum1$se3 + (e3 * parea)
  fsum1$se4 <- fsum1$se4 + e4
  fsum1$se5 <- fsum1$se5 + e5
  
  #CHECK THAT ADIMC >= UZTWC
  if (sma$adimc < sma$uztwc) {
    sma$adimc <- sma$uztwc
  }
  rsum <- rep(0, 7)
  #ADD TO SUMS OF RUNOFF COMPONENTS
  rsum[1] <- rsum[1] + sma$tlci
  rsum[2] <- rsum[2] + roimp
  rsum[3] <- rsum[3] + sdro
  rsum[4] <- rsum[4] + ssur
  rsum[5] <- rsum[5] + sif
  rsum[6] <- rsum[6] + bfs
  rsum[7] <- rsum[7] + bfp

  return(list(sma = sma, fsum1 = fsum1))
} # END OF FUNCTION FLAND1
