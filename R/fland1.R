fland1 <- function(sma, fsum1)
{
  EPS <- 1e-5
  zp <- sma$zperc
  
  edmnd <- sma$ep
  e1 <- edmnd * sma$uztwc / sma$uztwm
  red <- edmnd - e1
  sma$uztwc = sma$uztwc - e1
  if (abs(sma$uztwc) < EPS)
    sma$uztwc = 0
  e2 = 0
  if (sma$uztwc < 0) {
    e1 = e1 + sma$uztwc
    sma$uztwc = 0
    red = edmnd - e1
    if (sma$uzfwc < red) {
      e2 = sma$uzfwc
      sma$uzfwc = 0
      red <- red - e2
      # goto L225
    }
    else
    {
      e2 = red
      sma$uzfwc = sma$uzfwc - e2
      red = 0
    }
  }

  # Run in all cases except the goto statement above
  if(sma$uzfwc >= red){
      if ((sma$uztwc / sma$uztwm) < (sma$uzfwc / sma$uzfwm))
      {
        uzrat = (sma$uztwc + sma$uzfwc) / (sma$uztwm + sma$uzfwm)
        sma$uztwc = sma$uztwm * uzrat
        sma$uzfwc = sma$uzfwm * uzrat
      }
  }
  
  #L225
  e3 = red * sma$lztwc / (sma$uztwm + sma$lztwm)
  sma$lztwc = sma$lztwc - e3
  if (abs(sma$lztwc) < EPS)
    sma$lztwc = 0
  if (sma$lztwc < 0) {
    e3 = e3 + sma$lztwc
    sma$lztwc = 0
  }
  ratlzt = sma$lztwc / sma$lztwm
  saved = sma$rserv * (sma$lzfpm + sma$lzfsm)
  ratlz = (sma$lztwc + sma$lzfpc + sma$lzfsc - saved) / (sma$lztwm + sma$lzfpm +
                                                           sma$lzfsm - saved)
  
  if (ratlzt < ratlz) {
    del = (ratlz - ratlzt) * sma$lztwm
    
    sma$lztwc = sma$lztwc + del
    
    sma$lzfsc = sma$lzfsc - del
    if (sma$lzfsc < 0) {
      sma$lzfpc = sma$lzpfpc + sma$lzfsc
      sma$lzfsc = 0
    }
  }
  
  e5 = e1 + (red + e2) * (sma$adimc - e1 - sma$uztwc) / (sma$uztwm + sma$lztwm)
  
  sma$adimc = sma$adimc - e5
  if (abs(sma$adimc) < EPS)
    sma$adimc = 0
  if (sma$adimc < 0) {
    e5 = e5 + sma$adimc
    sma$adimc = 0
  }
  e5 = e5 * sma$adimp
  
  twx = sma$pxv + sma$uztwc - sma$uztwm
  if (twx < 0) {
    sma$uztwc = sma$uztwc + sma$pxv
    twx = 0
  }
  else {
    sma$uztwc = sma$uztwm
  }
  sma$adimc = sma$adimc + (sma$pxv - twx)
  
  roimp = sma$pxv * sma$pctim
  
  fsum1$simpvt = fsum1$simpvt - roimp
  
  sbf = ssur = sif = sperc = sdro = spbf = 0
  
  ninc = as.integer(1.0 + 0.2 * (sma$uzfwc + twx))
  if (ninc < sma$min_ninc)
    ninc = sma$min_ninc
  dinc = 1.0 /  ninc * sma$dt
  pinc = twx / ninc
  
  
  
  if (sma$uzk > 1)
    print(sma$uzk)
  if (sma$lzpk > 1)
    print(sma$lzpk)
  if (sma$lzsk > 1)
    print(sma$lzsk)
  duz = 1.0 - (1.0 - sma$uzk) ^ dinc
  dlzp = 1.0 - (1.0 - sma$lzpk) ^ dinc
  dlzs = 1.0 - (1.0 - sma$lzsk) ^ dinc
  parea = 1.0 - sma$adimp - sma$pctim
  
  for (i in 0:ninc) {
    adsur = 0
    ratio = (sma$adimc - sma$uztwc) / sma$lztwm
    addro = pinc * (ratio * ratio)
    sdro = sdro + (addro * sma$adimp)
    
    bf = sma$lzfpc * dlzp
    sma$lzfpc = sma$lzfpc - bf
    if (sma$lzfpc <= 1e-4) {
      bf = bf + sma$lzfpc
      sma$lzfpc = 0
    }
    sbf = sbf + bf
    spbf = spbf + bf
    bf = sma$lzfsc * dlzs
    sma$lzfsc = sma$lzfsc - bf
    if (sma$lzfsc <= 1e-4) {
      bf = bf + sma$lzfsc
      sma$lzfsc = 0
    }
    sbf = sbf + bf
    
    if ((pinc + sma$uzfwc) <= 1e-2) {
      sma$uzfwc = sma$uzfwc + pinc
      sma$adimc = sma$adimc + (pinc - addro - adsur)  #L249
      if (sma$adimc > (sma$uztwm + sma$lztwm))
      {
        addro = addro + (sma$adimc - (sma$uztwm + sma$lztwm))
        sma$adimc = sma$uztwm + sma$lztwm
      }
      next
    }
    percm = sma$lzfpm * dlzp + sma$lzfsm * dlzs
    
    if (zp < 0.0)
      zp = 0.0
    perc = percm * sma$uzfwc / sma$uzfwm
    
    defr = 1.0 - (sma$lztwc + sma$lzfpc + sma$lzfsc) / (sma$lztwm + sma$lzfpm + sma$lzfsm)
    
    if (defr < 0) {
      print(defr)
      c <- c(sma$lztwc, sma$lzfpc, sma$lzfsc)
      print(c)
      m <- c(sma$lztwm, sma$lzfpm, sma$lzfsm)
      print(m)
      stop("1")										#doubt
    }
    
    perc = perc * (1.0 + zp * (defr ^ sma$rexp))
    
    if (perc >= sma$uzfwc) {
      perc = sma$uzfwc
    }
    sma$uzfwc = sma$uzfwc - perc
    
    check = sma$lztwc + sma$lzfpc + sma$lzfsc + perc - sma$lztwm - sma$lzfpm -
      sma$lzfsm
    if (check > 0) {
      perc = perc - check
      sma$uzfwc = sma$uzfwc + check
    }
    
    sperc = sperc + perc
    del = sma$uzfwc * duz
    sif = sif + del
    sma$uzfwc = sma$uzfwc - del
    
    perct = perc * (1.0 - sma$pfree)
    if ((perct + sma$lztwc) <= sma$lztwm) {
      sma$lztwc = sma$lztwc + perct
      percf = 0
    }
    else {
      percf = perct + sma$lztwc - sma$lztwm
      sma$lztwc = sma$lztwm
    }
    
    percf = percf + (perc * sma$pfree)
    if (percf != 0) {
      hpl = sma$lzfpm / (sma$lzfpm + sma$lzfsm)
      
      ratlp = sma$lzfpc / sma$lzfpm
      ratls = sma$lzfsc / sma$lzfsm
      
      fracp = hpl * 2.0 * (1.0 - ratlp) / (1.0 - ratlp + 1.0 - ratls)
      if (fracp > 1.0)
        fracp = 1.0
      
      percp = percf * fracp
      percs = percf - percp
      sma$lzfsc = sma$lzfsc + percs
      if (sma$lzfsc > sma$lzfsm) {
        percs = percs + (-sma$lzfsc + sma$lzfsm)
        sma$lzfsc = sma$lzfsm
      }
      sma$lzfpc = sma$lzfpc + (percf - percs)
      
      if (sma$lzfpc > sma$lzfpm) {
        excess = sma$lzfpc - sma$lzfpm
        sma$lztwc = sma$lztwc + excess
        sma$lzfpc = sma$lzfpm
      }
    }
    
    if (pinc != 0) {
      if ((pinc + sma$uzfwc) <= sma$uzfwm) {
        sma$uzfwc = sma$uzfwc + pinc
      }
      else {
        sur = pinc + sma$uzfwc - sma$uzfwm
        sma$uzfwc = sma$uzfwm
        ssur = ssur + (sur * parea)
        adsur = sur * (1.0 - addro / pinc)
        ssur = ssur + (adsur * sma$adimp)
      }
    }
    sma$adimc = sma$adimc + (pinc - addro - adsur)  #L249
    if (sma$adimc > (sma$uztwm + sma$lztwm))
    {
      addro = addro + (sma$adimc - (sma$uztwm + sma$lztwm))
      sma$adimc = sma$uztwm + sma$lztwm
    }
  }
  eused = e1 + e2 + e3
  sif = sif * parea
  tbf = sbf * parea
  bfcc = tbf * (1.0 / (1.0 + sma$side))
  bfp = (spbf * parea) / (1.0 + sma$side)
  bfs = bfcc - bfp
  if (bfs < 0)
    bfs = 0
  bfncc = tbf - bfcc
  fsum1$sintft = fsum1$sintft + sif
  fsum1$sgwfp = fsum1$sgwfp + bfp
  fsum1$sgwfs = fsum1$sgwfs + bfs
  fsum1$srecht = fsum1$srecht + bfncc
  fsum1$srost = fsum1$srost + ssur
  fsum1$srodt = fsum1$srodt + sdro
  
  sma$tlci_flows[1] = roimp
  sma$tlci_flows[2] = sdro
  sma$tlci_flows[3] = ssur
  sma$tlci_flows[4] = sif
  sma$tlci_flows[5] = bfp
  sma$tlci_flows[6] = bfs
  sma$tlci_flows[7] = bfcc
  
  sma$tlci = roimp + sdro + ssur + sif + bfcc
  
  e4 = (edmnd - eused) * sma$riva
  
  sma$tlci = sma$tlci - e4
  
  if (sma$tlci < 0) {
    e4 = e4 + sma$tlci
    sma$tlci = 0
  }
  fsum1$srot = fsum1$srot + sma$tlci
  
  eused = eused * parea
  tet = eused + e5 + e4
  fsum1$sett = fsum1$sett + tet
  fsum1$se1 = fsum1$se1 + (e1 * parea)
  fsum1$se3 = fsum1$se3 + (e3 * parea)
  fsum1$se4 = fsum1$se4 + e4
  fsum1$se5 = fsum1$se5 + e5
  
  if (sma$adimc < sma$uztwc)
    sma$adimc = sma$uztwc
  rsum = rep(0, 7)
  rsum[1] = rsum[1] + sma$tlci
  rsum[2] = rsum[2] + roimp
  rsum[3] = rsum[3] + sdro
  rsum[4] = rsum[4] + ssur
  rsum[5] = rsum[5] + sif
  rsum[6] = rsum[6] + bfs
  rsum[7] = rsum[7] + bfp
  
  return(list(sma = sma, fsum1 = fsum1))
} #END OF FUNCTION FLAND1
