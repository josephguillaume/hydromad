update.filehashRDS<-function(object,...){
  switch(hydromad.getOption("parallel")[["update.runlist"]], 
         clusterApply = {
           parLapply(cl,names(object),function(n,object,...) {object[[n]]<-update(object[[n]],...);n},object=object,...)
         },
         lapply(names(object),function(n){object[[n]]<-update(object[[n]],...);n})
         )
  return(object) 
}
