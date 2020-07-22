
## Transfer function specifications of various orders
n0m0 <- c(v_s = 1.1)
n1m0 <- c(tau_s = 2)
n1m1 <- c(tau_s = 2, v_s = 0.9)
n2m0 <- c(tau_s = 30, tau_q = 2, series = 1) ## v_q = 1
n2m1 <- c(tau_s = 30, tau_q = 2, v_s = 0.3)
n2m2 <- c(tau_s = 30, tau_q = 2, v_s = 0.3, v_3 = 0.1)

## Third order models:
## "series = 0": three components in parallel (default)
n3s0 <- c(tau_s = 30, tau_q = 5, tau_3 = 1, v_s = 0.3, v_3 = 0.1) ## v_q = 1 - v_s - v_3
## "series = 1": two components in series and one in parallel
## (q & 3 are in series; s in parallel)
n3s1 <- c(tau_s = 30, tau_q = 5, tau_3 = 1, v_s = 0.3, v_3 = 0.7, series = 1) ## v_q = 1
## "series = 2": one component in series with two in parallel
## (3 in series; s & q in parallel)
n3s2 <- c(tau_s = 30, tau_q = 5, tau_3 = 1, v_s = 0.3, v_3 = 0.1, series = 2) ## v_q = 1 - v_s
## "series = 3": three components in series
n3s3 <- c(tau_s = 30, tau_q = 5, tau_3 = 1, v_3 = 0.7, series = 3) ## v_q = 1
