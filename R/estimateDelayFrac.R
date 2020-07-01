estimateDelayFrac <- function(DATA, rises = TRUE, lag.max = hydromad.getOption("max.delay")) {
  DATA <- as.ts(DATA)
  if (NROW(DATA) <= 1) {
    return(NA_integer_)
  }
  iQ <- 1
  if ("Q" %in% colnames(DATA)) {
    iQ <- which("Q" == colnames(DATA))[1]
  }
  Q <- DATA[, iQ]
  U <- DATA[, if (iQ == 1) {
    2
  } else {
    1
  }]
  if (all(Q == 0, na.rm = TRUE)) {
    return(NA_integer_)
  }
  do.rises <- rises
  if (do.rises) {
    rises <- function(x) {
      x[] <- c(0, pmax(diff(x), 0))
      x
    }
    Q <- rises(Q)
  }
  return(optimise(function(L) cor(lagFrac(U, L), Q, use = "complete.obs"), interval = c(0, lag.max), maximum = T)$maximum)
}

#################################################################################
# x<-mod$data
# event <- eventseq(x$P, thresh = 10, inthresh = 1, indur = 7, continue = TRUE)
# event.1lag <- lag(event, 1) # to consider the increment of the first obsQ
# whole.lag<-estimateDelayFrac(x)
#
# PQEM=mod$data
# lag.rain<-eventapply(PQEM,event.1lag,by.column=FALSE,function(x){
# ans.lag<-estimateDelayFrac(x, lag.max = 4)
# if (is.na(ans.lag)) {ans.lag<-whole.lag} # if obsQ=0 then delay=NA, put whole delay value
# x2<-merge(x, P1 = lagFrac(x[,1], ans.lag), ans.lag, all = c(TRUE, FALSE))
# return(x2)
# })

# require(reshape)
# melt.lag<-melt(lag.rain)
# write.csv(melt.lag, file = "lagrain.csv")
