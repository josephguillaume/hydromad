library(testthat)
library(hydromad)
library(patrick)

with_parameters_test_that(
  "zoo and hydromad handle date format: ",
  {
    if (!is.null(will_warn) && will_warn) {
      wrapper <- suppressWarnings
    } else {
      wrapper <- identity
    }
    
    DATA <- data.frame(
      P = c(100, rep(0, 9)),
      E = 20
    )

    DATA <- zoo(DATA, order.by = date)
    wrapper(expect_equal(nrow(as.ts(DATA)), 10))

    mod <- wrapper(hydromad(
      DATA,
      sma = "cwi",
      tw = 32, f = 2, scale = 0.01,
      routing = "expuh",
      tau_s = 1,
      warmup = 0
    ))
    Q <- c(
      63.2120558828558, 23.254415793483, 8.55482148687488,
      3.14714294791298, 1.15776918896487, 0.425919482241911, 0.156687021111184,
      0.0576419337652005, 0.0212052823815832, 0.00780098743241947
    )
    expect_equal(coredata(fitted(mod)), Q)
  },
  cases(
    `ymd` = list(date = sprintf("2020-01-%d", 1:10), will_warn = TRUE),
    `Date` = list(date = as.Date(sprintf("2020-01-%d", 1:10))),
    `POSIXct` = list(
      date = as.POSIXct(sprintf("2020-01-%d 3:00", 1:10), tz = "GMT"),
      will_warn = TRUE
    ),
    `POSIXct hourly` = list(
      date = as.POSIXct(sprintf("2020-01-01 %0d:00", 1:10), tz = "GMT"),
      will_warn = TRUE
    ),
    `POSIXlt` = list(date = as.POSIXlt(sprintf("2020-01-%d 3:00", 1:10), tz = "GMT"),will_warn = TRUE),
    `chron` = list(date = chron::as.chron(sprintf("2020-01-%d 3:00", 1:10)))
  )
)
