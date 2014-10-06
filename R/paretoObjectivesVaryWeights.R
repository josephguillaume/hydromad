paretoObjectivesVaryWeights<-function (MODEL, objective = hydromad.getOption("objective"), weights, fitBy,...)
{
  objective <- buildCachedObjectiveFun(objective, MODEL)
  switch(hydromad.getOption("parallel")[["paretoObjectivesVaryWeights"]],
         "clusterApply"={
           library(parallel)
           clusterExport(cl,c("objective","fitBy","MODEL"),envir=environment())
           front <- parApply(cl,weights,1,
                             ## fit model using weighted sum of objectives
                             function(w) fitBy(MODEL,objective=function(...) sum(w*sapply(objective,function(obj) obj(...))),
                                               ...))
         },
         front <- apply(weights,1,function(w) fitBy(MODEL,objective=function(...) sum(w*sapply(objective,function(obj) obj(...))),
                                                    ...))
         ) ##switch parallel
  names(front) <- apply(weights,1,paste,collapse="_")
  front <- as.runlist(front)
  return(front)
}
