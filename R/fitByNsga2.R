fitByNsga2<-function (MODEL, objective = hydromad.getOption("objective"), 
    control = hydromad.getOption("nsga2.control")) 
{
    start_time <- proc.time()
    objective <- buildCachedObjectiveFun(objective, MODEL)
    parlist <- as.list(coef(MODEL, warn = FALSE))
    isok <- sapply(parlist, function(x) !any(is.na(x)))
    parlist <- parlist[isok]
    isfixed <- (sapply(parlist, length) == 1)
    if (all(isfixed)) {
        warning("all parameters are fixed, so can not fit")
        return(MODEL)
    }
    thisVal <- objFunVal(thisMod, objective = objective)
    if (isTRUE(thisVal > bestFunVal)) {
      bestModel <<- thisMod
      bestFunVal <<- thisVal
    }
    return(-thisVal)
  }
  args <- modifyList(control, list(
    fn = do_nsga2, idim = length(parlist), odim = 1,
    lower.bounds = lower, upper.bounds = upper
  ))
  ans <- do.call(mco::nsga2, args)
  bestModel$funevals <- NA ## TODO
  bestModel$timing <- signif(proc.time() - start_time, 4)[1:3]
  bestModel$objective <- objective
  bestModel$fit.call <- match.call()
  bestModel$fit.result <- ans
  return(bestModel)
}
