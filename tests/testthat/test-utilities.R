library(testthat)
library(hydromad)

context("Utility functions")

test_that("flow conversions work correctly on known cases", {
  expect_equal(convertFlow(10, from = "cm", to = "in"), 3.93700787)
  expect_equal(convertFlow(0:10, from = "mm/day", to = "mm/hour"), (0:10) / 24)
  expect_equal(convertFlow(666, from = "m3/sec", to = "ML/day"), 666 * 24 * 60 * 60 / 1000)
  expect_equal(convertFlow(4, from = "mm", to = "ML", area.km2 = 2), 8)
  expect_error(convertFlow(1, from = "mm", to = "ML"))
  expect_equal(convertFlow(1, from = "ML / 15 minutes", to = "GL / year"), 4 * 24 * 365.25 / 1000)
})

test_that("parameterSets works correctly on different types of inputs", {
  set.seed(10)
  for (method in c("latin.hypercube", "random", "all.combinations")) {
    expect_equal(parameterSets(list(a = 10), 1, method = method)[, "a"], 10)
    expect_equal(parameterSets(list(a = c(10, 10)), 1, method = method)[, "a"], 10)
  }
  expect_equal(parameterSets(list(a = I(c(9, 10))), 1, method = "latin.hypercube")[, "a"], 9.5)
  expect_equal(parameterSets(list(a = I(c(9, 10))), 1, method = "random")[, "a"], 9)
})
