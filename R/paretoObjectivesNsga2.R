paretoObjectivesNsga2<-function (MODEL, objective = hydromad.getOption("objective"), 
    control = hydromad.getOption("nsga2.control")) 
{
    library(mco)
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
    parlist <- parlist[!isfixed]
    lower <- sapply(parlist, min)
    upper <- sapply(parlist, max)
    do_nsga2 <- function(pars) {
        names(pars) <- names(parlist)
        thisMod <- update(MODEL, newpars = pars)
        if (!isValidModel(thisMod)) 
            return(1e+08)
        thisVal <- objFunVal(thisMod, objective = objective)
        return(-unlist(thisVal))
    }
    args<-modifyList(control,list(
                                  fn=do_nsga2,
                                  idim=length(parlist),
                                  odim=length(objective),
                                  lower.bounds = lower, upper.bounds = upper
                                  ))
    ans <- do.call("nsga2",args)
    colnames(ans$par) <- names(parlist)
    MODEL$funevals <- NA ## TODO
    MODEL$timing <- signif(proc.time() - start_time, 4)[1:3]
    MODEL$objective <- objective
    MODEL$fit.call <- match.call()
    MODEL$fit.result <- ans
    if(isTRUE(hydromad.getOption("trace"))) cat("Creating runlist (running models on pareto front)\n")
    ## TODO: potential for parallelisation
    front <- as.runlist(apply(ans$par,1,function(p) update(MODEL,newpars=p)))
    return(front)
}
