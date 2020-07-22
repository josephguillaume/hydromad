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


#' Time-delay corrected performance measure
#'
#' adjVarTd coalesces modelled flow peaks to observed flow peaks for each event
#' separately. nseVarTd calculates Nash-Sutcliffe efficiency on the result
#' using \code{\link{nseStat}}.  Depending on the quality of the coalescing,
#' this better indicate performance ignoring timing error.
#'
#' The success of this method in minimising the effect of timing error depends
#' on how well modelled and observed peaks can be coalesced. This depends on:
#' \itemize{ \item \code{event} - The separation into events - too short events
#' result in spurious cross-correlations, too long events may not adequately
#' capture the variability in lag. Other settings of \code{\link{eventseq}} may
#' also have an effect. \item \code{lag.max} - How long a lag is considered.
#' Too long may result in correlations between peaks, too short will fail to
#' consider the true peak. Instead of passing it as an argument, consider
#' setting max.delay using \code{link{hydromad.options}} \item Other settings
#' of \code{\link{estimateDelay}} may also have an effect. \item The function
#' currently considers both positive and negative lag up to \code{lag.max}.
#' This can not be overridden. } Also note that large numbers of events will
#' run slower.
#'
#' @importFrom zoo zoo is.zoo index
#' @importFrom stats lag
#'
#' @aliases nseVarTd adjVarTd
#' @param obs observed data vector
#' @param mod model-predicted data vector corresponding to \code{obs}.
#' @param event zoo object of events, as returned by \code{\link{eventseq}}
#' @param \dots Additional arguments to \code{\link{nseStat}} and
#' \code{\link{estimateDelay}}.
#' @return For nseVarTd, a single numeric value.  For adjVarTd, a zoo object
#' with the original modelled and observed data, the adjusted model output and
#' the lag estimated for each event.
#' @author Joseph Guillaume
#' @seealso
#' \code{\link{hydromad.stats}},\code{\link{nseStat}},\code{\link{objFunVal}}
#' @keywords ts
#' @examples
#'
#'
#' data(Murrindindi)
#' x <- Murrindindi[1:100]
#' x <- merge(x, X = lag(x$Q, 2))
#'
#' event <- eventseq(x$P, thresh = 5, inthresh = 3.5, indur = 7, continue = TRUE)
#'
#' nseStat(x$Q, x$X)
#' nseVarTd(x$Q, x$X, event, lag.max = 3)
#'
#' ## Avoiding passing lag.max
#' hydromad.getOption("max.delay") ## Current setting - default is 10
#' hydromad.options(max.delay = 3)
#'
#' nseVarTd(x$Q, x$X, event)
#' hmadstat("r.sq.vartd")(x$Q, x$X, event = event)
#' @export
nseVarTd <- function(obs, mod, event, ...) {
  lag.mod <- adjVarTd(obs, mod, event, ...)
  return(nseStat(coredata(lag.mod$Q), coredata(lag.mod$mod1), ...))
}


#' @export
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
