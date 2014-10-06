paretoCatchments <- function(stat,show.models=NA,objectives="r.squared"){

  cat("\nCross-validation Pareto analysis\nWhich models cannot be rejected, due to dataset uncertainty/non-stationarity?\n")
  cat("\n== Eliminating Pareto-dominated models ==\n (worse than another model in all periods)\n")
  stat <- as.data.frame(stat)
  if(is.null(stat$Catchment)) stat$Catchment <- "Somewhere"
  stat.split <- split(stat,stat$Catchment)
  stat.split <- lapply(stat.split,areModelsDominated,objectives=objectives)

  
  nondom.long <- do.call("rbind",lapply(stat.split,getDominatedLong))
  stopifnot(!any(is.na(nondom.long$dominated)))
  nondom.long <- nondom.long[!nondom.long$dominated,]

  cat("\nHow many models are dominated (and therefore eliminated)\n")
  print(sapply(stat.split,function(res) table(res$dominated)))

  cat("\nWhich model structures are non-dominated in each catchment?\nProportion of model instances that are non-dominated\n")
  ##model.ok <- do.call(rbind,lapply(stat.split,function(xx) sort(unique(stat$Model.str)) %in% xx$Model.str[!xx$dominated]))
  model.ok <- do.call(rbind,lapply(stat.split,function(xx) {
    xx$Model.str <- factor(xx$Model.str)
    signif(table(xx$Model.str[!xx$dominated])/table(xx$Model.str),3)
  }))
  print(model.ok)
  
  if (!is.na(show.models)){
    isd <- do.call(rbind,stat.split)
    ## FIXME isd$calib.period <- period.labels[as.character(isd$calib.period)]
    isd$dominated <- ifelse(isd$dominated,"Yes","No")
    isd <- isd[with(isd,order(Catchment,dominated,Model.str,calib.period)),
               c(c("Catchment","dominated","Model.str","calib.period"),
                 setdiff(names(isd),c("Catchment","dominated","Model.str","calib.period")))
               ]
    if(isTRUE(show.models)){
      print(isd)
    } else if(!is.na(show.models)){
      cat("\nWriting csv of model realisations and whether dominated for each catchment\n")
      write.csv(isd,sprintf("%s_isdominated_models_catchments.csv",show.models),row.names=FALSE)
    }
  } else {
    cat("\nSpecify show.models=TRUE to show non-dominated and dominated models\nSpecify show.models=prefix to obtain csv of whether models are dominated\n")
  }
    
  cat("
== Performance across all periods ==
What is the range of non-dominated performance (RNDP) across all periods?
Is it large - is changing datasets causing problems?
")
  min.perf.all <- merge(
                        aggregate(value~Catchment,data=nondom.long,min),
                        aggregate(value~Catchment,data=nondom.long,max),
                        by=c("Catchment"),suffixes=c(".min",".max")
                        )
  min.perf.all$RNDP <- min.perf.all$value.max-min.perf.all$value.min
  print(min.perf.all)

    
  cat("
== Performance in each period ==
What is the RNDP in each period?
Is it low even though total RNDP is high? Why?
 Is there reason to believe the objective function is not comparable over time?
")
  min.perf.periods <- merge(
                            aggregate(value~sim.period+Catchment,data=nondom.long,min),
                            aggregate(value~sim.period+Catchment,data=nondom.long,max),
                            by=c("Catchment","sim.period"),suffixes=c(".min",".max")
                            )
  ## FIXME min.perf.periods$sim.period <- period.labels[as.character(min.perf.periods$sim.period)]
  min.perf.periods <- min.perf.periods[with(min.perf.periods,
                                            order(Catchment,sim.period)),]
  min.perf.periods$RNDP <- min.perf.periods$value.max-min.perf.periods$value.min
  print(min.perf.periods)

  cat("
== Worst non-dominated models in each period ==
Do any non-dominated models have unacceptable performance?
Which non-dominated model has the worst performance in each period? Why?
 Is it consistently the same dataset? Is there reason for that dataset to be problematic?
 Is it consistently the same model structure? Should another model structure have been selected?
 Is it consistently the same calibration objective function? Is it overfitting part of the record?
")
  for(cc in names(stat.split)){
    cat("\n",cc,"\n")
    wmodel <- apply(stat.split[[cc]][!stat.split[[cc]]$dominated,unique(stat$sim.period)],2,which.min)
    worst <- data.frame(sim.period=unique(stat$sim.period),
                        worst.performance=apply(stat.split[[cc]][!stat.split[[cc]]$dominated,unique(stat$sim.period)],2,min),
                        stat.split[[cc]][!stat.split[[cc]]$dominated,][wmodel,]
                        )
    worst$Catchment <- NULL         #redundant column - used as title
    worst$dominated <- NULL         #redundant column - subsetted to be false
    worst$objective <- NULL         #redundant column - specified as argument
    print(worst)
  }

  invisible(NULL)
}

getDominatedLong <- function(res.dom){
  resm2 <- as.data.frame(res.dom)
  resm2 <- melt(resm2,
                id.vars=intersect(names(res.dom),c("Model.str", "Catchment", "calib.period", "Cal.objfn","objective","dominated")),
                variable_name="sim.period")
  ## FIXME: ordering of sim.periods
  ##resm2$sim.period <- ordered(resm2$sim.period,levels=c("70","80","90","00"))
  ## FIXME: naming scheme of model realisations
  ##resm2$model <- paste(resm2$Model.str,resm2$calib.period,resm2$Cal.objfn)
  resm2
}
