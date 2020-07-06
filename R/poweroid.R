## poweroid geometry
## (c) Felix Andrews <felix@nfrac.org> 2009-06-04

## poweroid(V, alpha, beta) --> V, A, H, r
## poweroid(A, alpha, beta) --> V, A, H, r
## poweroid(H, alpha, beta) --> V, A, H, r
## poweroid(r, alpha, beta) --> V, A, H, r
## poweroid(V, A, alpha) --> beta
## poweroid(V, A, beta) --> alpha
## poweroid(V, r, alpha) --> beta
## poweroid(V, r, beta) --> alpha
## poweroid(V[:], A[:]) --> alpha, beta (fit)

poweroid <-
  function(alpha = NULL, beta = NULL,
           V = NULL, H = NULL,
           A = if (!missing(r)) pi * r^2,
           r = if (!missing(A)) sqrt(A / pi),
           ...,
           rel.error = FALSE,
           polish = FALSE, details = FALSE) {
    force(list(A, r))
    objective <-
      if (rel.error) {
        ~ mean(abs((obs - mod) / obs))
      } else {
        ~ mean(abs((obs - mod)))
      }
    if (!is.null(alpha) && !is.null(beta)) {
      ## geometry specified parametrically
      if (!is.null(V)) {
        A <- ((2 / beta + 1) * V / (alpha * pi^(-beta / 2)))^(2 / (2 + beta))
      } else if (!is.null(H)) {
        ## NOTE: this form may blow up with beta < ~0.01
        r <- (H / alpha)^(1 / beta)
        A <- pi * r^2
      }
      if (is.null(A)) {
        stop("give one of V, H, A or r")
      }
      ## now we have A, derive the rest
      r <- sqrt(A / pi)
      V <- (alpha * pi^(-beta / 2) / (2 / beta + 1)) * A^(1 + beta / 2)
      H <- alpha * r^beta
      return(data.frame(V = V, A = A, H = H, r = r))
    }
    if (!is.null(beta)) {
      ## beta (shape type) specified, derive alpha exactly
      if (!is.null(H) && !is.null(r)) {
        alpha <- H / (r^beta)
      } else
      if (!is.null(V) && !is.null(r)) {
        alpha <- (2 / beta + 1) * V / (pi * r^(2 + beta))
      } else
      if (!is.null(V) && !is.null(H)) {
        alpha <- (pi / ((2 / beta + 1) * V))^(beta / 2) * H^(1 + beta / 2)
        ## should avoid H ^ 2/beta
        ## alpha <- ((pi * H ^ (2/beta + 1)) / ((2 / beta + 1) * V)) ^ (beta/2)
      } else {
        stop("give any two of V, H, A/r to derive alpha")
      }
      return(list(alpha = alpha, beta = beta))
    }

    ## shape not specified; fit it to given data

    estAlpha <- function(betaEst) {
      ## find alpha values that fit data exactly, given this beta
      ## then take a (weighted) average
      alphaVec <- poweroid(V = V, H = H, A = A, beta = betaEst)$alpha
      if (rel.error) {
        alpha <- mean(alphaVec)
      } else {
        w <- if (!is.null(A)) A else H
        alpha <- weighted.mean(alphaVec, w = w)
      }
      alpha
    }
    objfn <- function(x, withAlpha = FALSE) {
      if (withAlpha) {
        alpha <- x[1]
        beta <- x[2]
      } else {
        beta <- x[1]
        alpha <- estAlpha(beta)
      }
      if (!is.null(H) && !is.null(A)) {
        H.hat <- alpha * r^beta
        foo <- list(obs = H, mod = H.hat)
      } else
      if (!is.null(V) && !is.null(A)) {
        A.hat <- ((2 / beta + 1) * V / (alpha * pi^(-beta / 2)))^(2 / (2 + beta))
        foo <- list(obs = A, mod = A.hat)
      } else
      if (!is.null(V) && !is.null(H)) {
        H.hat <- ((2 / beta + 1) * V / pi)^(beta / (2 + beta)) * alpha^(2 / (2 + beta))
        foo <- list(obs = H, mod = H.hat)
      } else {
        stop("need more information")
      }
      objtmp <- objective
      if (inherits(objective, "formula")) {
        objtmp <- objective[[2]]
      }
      eval(objtmp, foo)
    }
    ## run optimize with different starting points for beta
    results <- lapply(
      c(0.5, 1, 1.5, 2, 2.5, 3.5),
      function(beta) {
        optimize(objfn, interval = pmax(0.01, beta + c(-0.5, 0.5)), ...)
      }
    )
    if (details) {
      return(results)
    }
    ## extract overall best result (value of beta)
    whichBest <- which.min(sapply(results, function(x) x$objective))
    beta <- results[[whichBest]]$minimum
    alpha <- estAlpha(beta)
    ## now polish result
    ## this does not seem to change the result at all:
    if (polish) {
      result <- nlminb(c(alpha = alpha, beta = beta), objfn,
        lower = c(0, 0.01),
        withAlpha = TRUE
      )
      if (result$convergence != 0) {
        warning(result$message)
      }
      alpha <- result$par[1]
      beta <- result$par[2]
    }
    return(list(alpha = alpha, beta = beta))
    #    result <- nlm(objfn, c(alpha = alpha, beta = beta),
    #                     withAlpha = TRUE, ...)
    #    if (result$code == 4) {
    #        warning("iteration limit exceeded")
    #    } else if (result$code > 3) {
    #        warning("code ", result$code, ", see ?nlm")
    #    }
    #    return(list(alpha = result$estimate[1],
    #                beta = result$estimate[2]))
  }
