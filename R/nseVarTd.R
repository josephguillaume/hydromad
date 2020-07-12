## NSE with variable time delay corrected Qmod

## Lag may be positive or negative

## nseVarNonintTd <- function(obs,mod,event,...){
##     event.1lag <- lag(event, 1) # to consider the increment of the first obsQ
##     whole.lag<-estimateDelayFrac(x)

##     PQEM=cbind(U=mod,Q=obs)
##     lag.rain<-eventapply(PQEM,event.1lag,by.column=FALSE,function(x){
##         ans.lag<-estimateDelayFrac(x, lag.max = 4)
##         if (is.na(ans.lag)) {ans.lag <- whole.lag} # if obsQ=0 then delay=NA, put whole delay value
##         x2<-merge(x, mod1 = lagFrac(PQEM[,1], ans.lag), all = c(TRUE, FALSE))
##         return(x2)
##     })
##     return(nseStat(coredata(lag.rain$Q), coredata(lag.rain$mod1), ...))
## }

adjVarTd <- function(obs, mod, event, ...) {
  ## Needs to be zoo objects with same indices
  stopifnot(is.zoo(event))
  stopifnot(identical(index(obs), index(mod)))
  stopifnot(identical(index(obs), index(event)))
  ## Shift events by one day to allow rises to be better picked up
  event.1lag <- merge(lag(event, 1), zoo(, index(event)))
  event.1lag[length(event.1lag)] <- event.1lag[length(event.1lag) - 1]

  mod.obs <- cbind(U = mod, Q = obs)
  whole.lag <- estimateDelay(mod.obs, negative.ok = T, ...)
  lag.mod <- eventapply(mod.obs, event.1lag, by.column = FALSE, function(x) {
    ans.lag <- estimateDelay(x, negative.ok = T, ...)
    if (is.na(ans.lag)) {
      ans.lag <- whole.lag
    } # if obsQ=0 then delay=NA, put whole delay value
    x2 <- merge(x, mod1 = lag(mod.obs[, 1], ans.lag), ans.lag, all = c(TRUE, FALSE))
    return(x2)
  })
  return(do.call(rbind, lag.mod))
}


nseVarTd <- function(obs, mod, event, ...) {
  lag.mod <- adjVarTd(obs, mod, event, ...)
  return(nseStat(coredata(lag.mod$Q), coredata(lag.mod$mod1), ...))
}

################################################################################

## library(hydromad)

## data(Murrindindi)
## x <- Murrindindi[1:100]
## x <- merge(x,X=lag(x$Q,2))

## event <- eventseq(x$P, thresh = 5, inthresh = 3.5, indur = 7, continue = TRUE)
## table(coredata(event))
## ## Will crash if there's too few events

## nseStat(x$Q,x$X)
## nseVarTd(x$Q,x$X,event)
## ##nseVarNonintTd(x$Q,x$X,event)
## hmadstat("r.sq.vartd")(x$Q,x$X,event=event)
