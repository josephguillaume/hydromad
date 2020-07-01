runRSM <- function(modx, ..., objective = hydromad.getOption("objective")) {
  if (!requireNamespace("rsm")) stop("package rsm is required for runRSM")
  ## Sample using design
  bounds <- getFreeParsRanges(modx)

  coding <- lapply(
    names(bounds),
    function(var) {
      max.var <- bounds[[var]][1]
      min.var <- bounds[[var]][2]
      as.formula(sprintf("%s.c ~ (2*%s-(%f+%f))/(%f-%f)", var, var, max.var, min.var, max.var, min.var))
    }
  )

  evals <- rsm::ccd(as.formula(sprintf("~%s", paste(sprintf("%s.c", names(bounds)), collapse = "+"))),
    coding = coding, inscribed = TRUE, ...
  )

  cat(sprintf("Running %d model evaluations\n", nrow(evals)))
  evals$Y <- evalPars(decode.data(evals)[, names(bounds)], modx, objective = objective)

  ################################################################################
  ## Estimate quadratic response surface
  form <- as.formula(sprintf("Y~SO(%s)", paste(sprintf("%s.c", names(bounds)), collapse = ",")))
  evals.rsm <- rsm::rsm(form, data = evals)

  return(evals.rsm)
}

evalRSM <- function(modx, rsm.object, n = 100,
                    objective = hydromad.getOption("objective"),
                    method = "latin.hypercube") {
  indep.sample <- parameterSets(getFreeParsRanges(modx), samples = n, method = method)
  indep.sample.coded <- coded.data(indep.sample, formulas = codings(rsm.object))
  cat(sprintf("Running %d model evaluations\n", nrow(indep.sample)))
  indep.sample.coded$Y <- evalPars(decode.data(indep.sample), modx, objective = objective)
  f <- predict(rsm.object, newdata = indep.sample.coded)
  r <- indep.sample.coded$Y - f
  ans <- list(
    fitted.values = f,
    model.values = indep.sample.coded$Y,
    indep.sample = indep.sample.coded
  )
  rss <- sum(r^2)
  mss <- sum((f - mean(f))^2)
  rdf <- n - rsm.object$rank
  ans$df <- c(rsm.object$rank, rdf, NCOL(rsm.object$qr$qr))
  resvar <- rss / rdf
  ans$sigma <- sqrt(resvar)
  ans$r.squared <- mss / (mss + rss)
  df.int <- 1 ## intercept present
  ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int) / rdf)
  class(ans) <- "summary.rsm.hydromad"
  ans
}

print.summary.rsm.hydromad <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  rdf <- x$df[2L]
  cat("\nResidual standard error:", format(signif(
    x$sigma,
    digits
  )), "on", rdf, "degrees of freedom")
  cat("\n")
  cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
  cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits))
  cat("\n")
}
