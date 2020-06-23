
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hydromad

<!-- badges: start -->

<!-- badges: end -->

The goal of hydromad is to provide a modelling framework for
environmental hydrology: water balance accounting and flow routing in
spatially aggregated catchments.

Hydromad supports simulation, estimation, assessment and visualisation
of flow response to time series of rainfall and other drivers. A minimal
unit hydrograph framework is used, where areal rainfall is passed
through a soil moisture accounting (SMA) model to estimate effective
rainfall; this is then passed through a routing model to estimate
streamflow. Included are several implementations of common hydrological
models consistent with this framework.

## Installation

<!--- You can install the released version of hydromad from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("hydromad")
```
--->

Currently you can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JosephGuillaume/hydromad")
```

## Usage

This is a basic example which shows you how to solve a common problem:

``` r
library(hydromad)
#> Loading required package: zoo
#> 
#> Attaching package: 'zoo'
#> The following objects are masked from 'package:base':
#> 
#>     as.Date, as.Date.numeric
#> Loading required package: lattice
#> Loading required package: latticeExtra
#> Loading required package: polynom
#> Loading required package: reshape
## basic example code
data(Cotter)
## IHACRES CWI model with exponential unit hydrograph
## an unfitted model, with ranges of possible parameter values
modx <- hydromad(Cotter[1:1000], sma = "cwi", routing = "expuh",
                 tau_s = c(2,100), v_s = c(0,1))
```
