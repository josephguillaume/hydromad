## Load hydromad package
library(hydromad)
## Load sensitivity package
library(sensitivity)

## Load data
data(Cotter)
obs <- Cotter[1:1000]

## Define rainfall-runoff model structure
modx <- hydromad(obs,
  sma = "cwi", routing = "expuh",
  tau_s = c(2, 100), v_s = c(0, 1)
)

## Set the random seed to obtain replicable results
set.seed(19)

################################################################################
## Sensitivity using Morris method of NSE* objective function to
##  IHACRES-CWI model parameters using a subset of data from Cotter catchment

## Run Morris Sensitivity analysis
mm <- morris(
  ## Hydromad's Function to evaluate a target function for varying parameters
  model = evalPars,
  ## names of parameters to be varied
  factors = names(getFreeParsRanges(modx)),
  ## number of elementary effects computed per factor
  r = 4,
  design = list(
    ## Uses Morris (1991) One-at-a-Time (OAT) sampling scheme (see ?morris)
    type = "oat",
    ## Number of levels of parameter values
    levels = 10,
    ## Number of levels that are increased/decreased for computing the elementary effects
    grid.jump = 2
  ),
  ## Minimum value for each parameter
  binf = sapply(getFreeParsRanges(modx), min),
  ## Maximum value for each parameter
  bsup = sapply(getFreeParsRanges(modx), max),
  ## Arguments to be passed to evalPars: model and objective
  object = modx,
  ## NSE* objective function
  objective = ~ hmadstat("r.squared")(Q, X) / (2 - hmadstat("r.squared")(Q, X))
)

print(mm)

## Default plot of mu.star (mean absolute sensitivity) and sigma (sd of sensitivity)
plot(mm, main = "Sensitivity NSE*~IHACRES-CWI parameters with Cotter data")

## For custom plots, mu.star and sigma can be explicitly calculated
mu.star <- apply(mm$ee, 2, function(x) mean(abs(x)))
sigma <- apply(mm$ee, 2, sd)
plot(mu.star, sigma,
  pch = 20, xlab = expression(mu^"*"),
  ylab = expression(sigma)
)
text(mu.star, sigma, labels = colnames(mm$ee), pos = 4, offset = 0.4)

################################################################################
## We then perform sensitivity analysis on the prediction function, to identify
## whether parameters that are insensitive with the objective function can be fixed:
##  Sensitivity using Morris method of F20 prediction function
##  (frequency/ratio of days of flow below 20%ile)
##  to IHACRES-CWI model parameters using a subset of data from Cotter catchment

## Calculate observed 20%ile flow
thres.Q20 <- as.numeric(quantile(obs$Q, probs = c(0.2), na.rm = TRUE))
## Define function to calculate number of days below thres.Q20
F20 <- function(X) length(which(X < thres.Q20)) / length(X)

mm <- morris(
  model = evalPars,
  factors = names(getFreeParsRanges(modx)),
  r = 4,
  design = list(
    type = "oat",
    levels = 10,
    grid.jump = 2
  ),
  binf = sapply(getFreeParsRanges(modx), min),
  bsup = sapply(getFreeParsRanges(modx), max),
  object = modx,
  ## Change objective to use the prediction function defined above
  ## See ?objFunVal
  ## The objective function can refer to Q and X,
  ##  representing observed and modelled flow, respectively.
  ##  It should return a single numeric value.
  objective = ~ F20(X)
)

print(mm)
plot(mm, main = "Sensitivity F20~IHACRES-CWI parameters with Cotter data")

## See the Morris (1991) journal paper(?morris) for a more comprehensive interpretation.
## These results are only an example, but if sampling were sufficient, the results would
##  suggest that NSE* is insensitive to tau_s (mu.star=0.065), but F20 is sensitive to
##  tau_s (mu.star=0.34)
## The parameter value selected therefore affects the prediction, so tau_s cannot be
##  fixed to improve identifiability unless its value is sufficiently certain.
## Using a different objective function, e.g. NSElog* might improve identifiability.
## The settings for Morris used may also influence results. Try changing: the parameter ranges,
##  number of replicates, number of levels and length of data series.
## A SOBOL sensitivity analysis would also provide confidence intervals around the TSI.


################################################################################
## Sensitivity using SOBOL2002 method of NSElog* prediction function to
##  IHACRES-CWI model parameters using a subset of data from Cotter catchment
## This might take a while, potentially ~15min

## Save current time
st <- proc.time()

## Set number of samples desired
n <- 1000
ss <- sobol2002(
  model = evalPars,
  ## Draw two random samples of parameters
  X1 = parameterSets(getFreeParsRanges(modx), n),
  X2 = parameterSets(getFreeParsRanges(modx), n),
  ## Number of bootstrap replicates
  nboot = 100,
  ## Arguments to be passed to evalPars
  object = modx,
  ## NSElog* objective function (using the logarithm of Q and X)
  objective = ~ hmadstat("r.sq.log")(Q, X) / (2 - hmadstat("r.sq.log")(Q, X))
)

## Print elapsed time
print(proc.time() - st)

## Show results
print(ss)
plot(ss)

## tau_s appears to still be insensitive with NSElog* (TSI=0.127), suggesting
## the information in this dataset may not be sufficient to identify this parameter.
