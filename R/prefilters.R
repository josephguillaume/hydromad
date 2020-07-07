## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##


defaultPrefilters <- function() {
  alpha_q <- c(0.01, 0.2)
  alpha_s <- c(0.9, 0.95, 0.98, 0.8)

  alpha <- expand.grid(alpha_q = alpha_q, alpha_s = alpha_s)

  a1 <- with(alpha, alpha_s + alpha_q)
  a2 <- with(alpha, -(alpha_s * alpha_q))
  aa <- data.frame(a_1 = a1, a_2 = a2)
  split(as.matrix(aa), 1:NROW(aa))
}

makePrefilter <-
  function(DATA,
           order = hydromad.getOption("order"),
           pureAR = FALSE,
           na.action = na.exclude) {
    ## get data into the right form
    if (!is.ts(DATA)) DATA <- as.ts(DATA)
    if (NCOL(DATA) > 1) {
      DATA <- DATA[, "Q"]
    }
    DATA <- na.action(DATA)
    ## check values
    stopifnot(length(order) == 2)
    order[1] <- max(1, order[1])
    ## generate prefilter from data
    ## pure AR method is very fast
    if (pureAR) {
      return(ar(DATA,
        order.max = order[1], demean = FALSE,
        na.action = na.action
      )$ar)
    }
    ## otherwise fit an ARMA model
    mod <- arima(DATA,
      order = c(order[1], 0, order[2] + 1),
      include.mean = FALSE
    )
    coef(mod)[1:(order[1])]
  }


alpha.init <-
  function(n = 2) {
    stopifnot(n %in% 1:2)

    alpha_q <- c(0.01, 0.2)
    alpha_s <- c(0.9, 0.95, 0.98, 0.8)

    if (n == 1) {
      return(data.frame(alpha_s = alpha_q))
    }

    alpha <- expand.grid(alpha_q = alpha_q, alpha_s = alpha_s)
    alpha <- subset(alpha, alpha_q <= alpha_s)
    return(alpha)
    # unlist(alpha[k,])
  }

a.init <-
  function(n = 2) {
    alpha <- alpha.init(n = n)
    if (n == 1) {
      return(alpha)
    }
    a1 <- with(alpha, alpha_s + alpha_q)
    a2 <- with(alpha, -(alpha_s * alpha_q))
    data.frame(a_1 = a1, a_2 = a2)
  }
