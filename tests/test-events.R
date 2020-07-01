library(testthat)
library(hydromad)

context("Event handling functions")

data(SalmonBrook)
dat <- window(SalmonBrook, start = "1990-01-01", end = "1994-01-01")

evp <- eventseq(dat$P, mingap = 7)
evq <- eventseq(dat$Q, thresh = 3)
evq.ts <- eventseq(as.ts(dat$Q), thresh = 3)
evpq <- eventseq(dat$P,
  thresh = 5, mindur = 3,
  inx = dat$Q, inthresh = 0.5
)

test_that("eventseq seems to work", {
  expect_that(evp, is_a("zoo"))
  expect_that(evq, is_a("zoo"))
  expect_that(evq.ts, is_a("zoo"))
  expect_that(evpq, is_a("zoo"))
  expect_that(coredata(evp), is_a("factor"))
  expect_that(coredata(evpq), is_a("factor"))
  expect_that(index(evp), equals(index(dat)))
  expect_that(index(evq.ts), equals(index(as.ts(dat))))
  expect_that(c(unclass(evq)), equals(c(unclass(evq.ts))))
  expect_that(nlevels(evp), equals(39))
  expect_that(nlevels(evq), equals(14))
  expect_that(sum(is.na(coredata(evp))), equals(524))
})

test_that("findThresh seems to work", {
  set.seed(0)
  x <- rnorm(100)
  t1 <- findThresh(x, n = 20)
  t2 <- findThresh(x, n = 5, mingap = 2)
  expect_that(abs(nlevels(eventseq(x, t1)) - 20) <= 2, is_true())
  expect_that(nlevels(eventseq(x, t2, mingap = 2)) - 5, equals(0))
})

test_that("eventapply seems to work with single series", {
  ## NOTE need to handle functions returning vectors as well as scalars.
  ## (1) scalar result:
  psums <- eventapply(dat$P, evp)
  expect_that(psums, is_a("zoo"))
  expect_that(index(psums), is_a("Date"))
  expect_that(NCOL(psums), equals(1))
  expect_that(NROW(psums), equals(nlevels(evp)))
  ## factor events are not sync'd (cbinded) with the data series
  ## but here we know that they are already synchronised.
  expect_that(
    as.vector(psums),
    is_identical_to(as.vector(eventapply(dat$P, coredata(evp))))
  )
  ## (2) vector (>1) result:
  p2num <- eventapply(dat$P, evp,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  expect_that(colnames(p2num), equals(c("mean", "sd")))
  expect_that(NROW(p2num), equals(nlevels(evp)))
  ## variable length result, with simplify = FALSE
  pvari <- eventapply(dat$P, evp, FUN = coredata, simplify = FALSE)
  expect_that(pvari, is_a("list"))
  expect_that(length(pvari), equals(nlevels(evp)))
  expect_that(names(pvari), equals(format(unname(index(psums)))))
})

test_that("eventapply seems to work with multiple series", {
  ## (1) scalar result with by.column = FALSE
  durs <- eventapply(dat, evp, FUN = nrow, by.column = FALSE)
  expect_that(durs, is_a("zoo"))
  expect_that(NCOL(durs), equals(1))
  expect_that(NROW(durs), equals(nlevels(evp)))
  ## (2) scalar result with by.column = TRUE (the default)
  sums <- eventapply(dat, evp, FUN = sum)
  expect_that(sums, is_a("zoo"))
  expect_that(NCOL(sums), equals(NCOL(dat)))
  expect_that(colnames(sums), equals(colnames(dat)))
  ## (3) vector result with by.column = FALSE
  ## should be exactly the same in this case!
  sums2 <- eventapply(dat, evp, FUN = colSums, by.column = FALSE)
  expect_that(sums, equals(sums2))
  ## (4) vector result with by.column = TRUE
  each2num <- eventapply(dat, evp,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  expect_that(
    colnames(each2num),
    equals(c("P.mean", "P.sd", "Q.mean", "Q.sd", "E.mean", "E.sd"))
  )
  expect_that(index(each2num), equals(index(sums)))
})

test_that("eventinfo looks ok", {
  info <- eventinfo(dat$P, evp)
  expect_that(info, is_a("data.frame"))
  expect_that(
    colnames(info),
    is_identical_to(c(
      "Time", "Month", "Year",
      "Value", "Duration", "PreDuration"
    ))
  )
  expect_that(nrow(info), equals(nlevels(evp)))
})
