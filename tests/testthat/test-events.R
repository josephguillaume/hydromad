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
  expect_is(evp, "zoo")
  expect_is(evq, "zoo")
  expect_is(evq.ts, "zoo")
  expect_is(evpq, "zoo")
  expect_is(coredata(evp), "factor")
  expect_is(coredata(evpq), "factor")
  expect_equal(index(evp), index(dat))
  expect_equal(index(evq.ts), index(as.ts(dat)))
  expect_equal(c(unclass(evq)), c(unclass(evq.ts)))
  expect_equal(nlevels(evp), 39)
  expect_equal(nlevels(evq), 14)
  expect_equal(sum(is.na(coredata(evp))), 524)
})

test_that("findThresh seems to work", {
  set.seed(0)
  x <- rnorm(100)
  t1 <- findThresh(x, n = 20)
  t2 <- findThresh(x, n = 5, mingap = 2)
  expect_lte(abs(nlevels(eventseq(x, t1)) - 20), 2)
  expect_equal(nlevels(eventseq(x, t2, mingap = 2)) - 5, 0)
})

test_that("eventapply seems to work with single series", {
  ## NOTE need to handle functions returning vectors as well as scalars.
  ## (1) scalar result:
  psums <- eventapply(dat$P, evp)
  expect_is(psums, "zoo")
  expect_is(index(psums), "Date")
  expect_equal(NCOL(psums), 1)
  expect_equal(NROW(psums), nlevels(evp))
  ## factor events are not sync'd (cbinded) with the data series
  ## but here we know that they are already synchronised.
  expect_identical(
    as.vector(psums),
    as.vector(eventapply(dat$P, coredata(evp)))
  )
  ## (2) vector (>1) result:
  p2num <- eventapply(dat$P, evp,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  expect_equal(colnames(p2num), c("mean", "sd"))
  expect_equal(NROW(p2num), nlevels(evp))
  ## variable length result, with simplify = FALSE
  pvari <- eventapply(dat$P, evp, FUN = coredata, simplify = FALSE)
  expect_output(str(class(pvari)), "list")
  expect_equal(length(pvari), nlevels(evp))
  expect_equal(names(pvari), format(unname(index(psums))))
})

test_that("eventapply seems to work with multiple series", {
  ## (1) scalar result with by.column = FALSE
  durs <- eventapply(dat, evp, FUN = nrow, by.column = FALSE)
  expect_is(durs, "zoo")
  expect_equal(NCOL(durs), 1)
  expect_equal(NROW(durs), nlevels(evp))
  ## (2) scalar result with by.column = TRUE (the default)
  sums <- eventapply(dat, evp, FUN = sum)
  expect_is(sums, "zoo")
  expect_equal(NCOL(sums), NCOL(dat))
  expect_equal(colnames(sums), colnames(dat))
  ## (3) vector result with by.column = FALSE
  ## should be exactly the same in this case!
  sums2 <- eventapply(dat, evp, FUN = colSums, by.column = FALSE)
  expect_equal(sums, sums2)
  ## (4) vector result with by.column = TRUE
  each2num <- eventapply(dat, evp,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  expect_equal(
    colnames(each2num),
    c("P.mean", "P.sd", "Q.mean", "Q.sd", "E.mean", "E.sd")
  )
  expect_equal(index(each2num), index(sums))
})

test_that("eventinfo looks ok", {
  info <- eventinfo(dat$P, evp)
  expect_output(str(class(info)), "data.frame")
  expect_identical(
    colnames(info),
    c(
      "Time", "Month", "Year",
      "Value", "Duration", "PreDuration"
    )
  )
  expect_equal(nrow(info), nlevels(evp))
})
