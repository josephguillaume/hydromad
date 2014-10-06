## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Joseph Guillaume <josephguillaume@gmail.com>
##

crossValidate <- function(MODEL,periods,
                          name.Model.str=paste(MODEL$sma,MODEL$routing),
                          name.Cal.objfn="unknown",
                          fitBy,...){
  cv.set <- splitData(MODEL,periods=periods)
  runs <- runlist()                   
  for(n in names(cv.set)){
    if(isTRUE(hydromad.getOption("trace"))) cat("\nFitting period: ",n,"\n")    
    fitx <- fitBy(cv.set[[n]],...)
    new.runs <- update(cv.set,newpars=coef(fitx))
    names(new.runs) <- sprintf("%s_cal%s",names(cv.set),n)
    ## Preserve fit attributes
    new.runs[[sprintf("%s_cal%s",n,n)]] <- fitx
    for(m in names(new.runs)){
      ##new.runs[[m]]$objective <- fitx$objective
      new.runs[[m]]$name.Model.str <- name.Model.str
      new.runs[[m]]$name.Cal.objfn <- name.Cal.objfn
    }
    runs <- c(runs,new.runs)
  }
  class(runs) <- unique(c("crossvalidation",class(runs),"runlist", "list"))
  runs
}

summary.crossvalidation <- function(object, ...){
  s=NextMethod(object,...)
  s$sim.period <- sub("_cal.*","",rownames(s))
  s$calib.period <- sub(".*_cal","",rownames(s))
  s$Model.str <- sapply(object,function(x) x$name.Model.str)
  s$Cal.objfn <- sapply(object,function(x) x$name.Cal.objfn)
  s
}  
