paretoTimeAnalysis <- function(rl,
                             objectives="r.squared",
                             qoi=list(
                               q90=function(Q,X,...) quantile(X,0.9,na.rm=TRUE),
                               r.sq.log=hmadstat("r.sq.log")
                               ),...
                             ){
  stopifnot(inherits(rl,"crossvalidation"))
  stat <- summary(rl,...)
  paretoCatchments(stat,objectives=objectives)

  cat("
== Performance with other statistics ==
Use the set of non-dominated models as an ensemble to predict a quantity of interest
Is the model unacceptable in any period? Is the uncertainty too large?
")

  ## TODO: avoid recomputation of dominated models cf paretoCatchments
  stat <- as.data.frame(stat)
  if(is.null(stat$Catchment)) stat$Catchment <- "Somewhere"
  stat.split <- split(stat,stat$Catchment)
  stat.split <- do.call(rbind,lapply(stat.split,areModelsDominated,objectives=objectives))
  id.vars <- intersect(names(stat), c("Model.str", "Catchment", 
                                     "calib.period", "Cal.objfn"))
  is.nondominated <- apply(stat[,id.vars],1,paste,collapse="_") %in% apply(stat.split[!stat.split$dominated,id.vars],1,paste,collapse="_")
  stopifnot(length(is.nondominated)==length(rl))
  
  ## Calculate qois for each model
  p <- lapply(rl[is.nondominated],objFunVal,objective=qoi)
  ## Convert from list of lists to data.frame
  p <- as.data.frame(do.call(rbind,lapply(p,unlist)))
  ## Extract simulation and calibration period names
  p$sim.period <- sub("_cal.*", "", rownames(p))
  p$calib.period <- sub(".*_cal", "", rownames(p))
  ## Melt qois
  pm <- melt(p,id.vars=intersect(names(p),c("Model.str", "Catchment", "sim.period","calib.period", "Cal.objfn","objective","dominated")))
  ## Cast and aggregate showing min, max and range of each qoi
  p2 <- cast(pm,variable+sim.period~.,
             fun.aggregate=function(x) c(min=min(x),max=max(x),range=diff(range(x))))
  print(p2)

  invisible(NULL)
}##paretoTimeAnalysis
