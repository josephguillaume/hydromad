findUnivariateBounds <- function(modx, fitx, thres, objective = hydromad.getOption("objective")) {
  pars <- getFreeParsRanges(modx)
  for (p in names(pars)) {
    best <- coef(fitx)[[p]]
    ## Lower side
    try({
      lo <- uniroot(function(val) {
        names(val) <- p
        objFunVal(update(fitx, newpars = val), objective = objective) - thres
      }, interval = c(min(pars[[p]]), best))
      pars[[p]][1] <- lo$root
    })
    ## upper side
    try({
      hi <- uniroot(function(val) {
        names(val) <- p
        objFunVal(update(fitx, newpars = val), objective = objective) - thres
      }, interval = c(best, max(pars[[p]])))
      pars[[p]][2] <- hi$root
    })
  }
  return(pars)
}
