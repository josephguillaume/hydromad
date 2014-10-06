## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Joseph Guillaume <josephguillaume@gmail.com>
##

splitData<-function(object,start.dates=NULL,end.dates=NULL,periods=NULL){
  ## Split data according to start.dates and end.dates
  if(!is.null(periods)){
    start.dates <- sapply(periods,head,1)
    end.dates <- sapply(periods,tail,1)
  }
  cv.data <- mapply(window,start=start.dates,end=end.dates,MoreArgs=list(x=observed(object,all=TRUE,select=TRUE)),SIMPLIFY=FALSE)
  rl=as.runlist(mapply(update,newdata=cv.data,MoreArgs=list(object=object),SIMPLIFY=FALSE))
  if(!is.null(names(periods))) names(rl)<-names(periods)
  else names(rl) <- paste(start.dates,end.dates,sep="_")
  return(rl)
}

