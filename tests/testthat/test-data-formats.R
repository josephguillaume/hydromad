library(testthat)
library(hydromad)
library(patrick)

with_parameters_test_that(
  "zoo and hydromad handle date format: ",
  {
    DATA <- data.frame(
      P = c(100, rep(0, 9)),
      E = 20
    )

    DATA <- zoo(DATA, order.by = date)
    # TODO: check for unexpected warnings?
    # In some cases, expect ‘x’ does not have an underlying regularity
    suppressWarnings(expect_equal(nrow(as.ts(DATA)), 10))

    for (w in hydromad_warning) {
      expect_warning(
        {
          mod <- hydromad(
            DATA,
            sma = "cwi",
            tw = 32, f = 2, scale = 0.01,
            routing = "expuh",
            tau_s = 1,
            warmup = 0
          )
        },
        w
      )
    }

    # TODO: check for unexpected warnings?
    suppressWarnings(mod <- hydromad(
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
    `ymd` = list(
      date = sprintf("2020-01-%d", 1:10),
      as.ts_warning = c("‘x’ does not have an underlying regularity"),
      hydromad_warning = c(
        "DATA appears to cover less than 60 days.*"
      )
    ),
    `Date` = list(date = as.Date(sprintf("2020-01-%d", 1:10)), hydromad_warning = c(
      "DATA appears to cover less than 60 days.*"
    )),
    `POSIXct` = list(
      date = as.POSIXct(sprintf("2020-01-%d 3:00", 1:10), tz = "GMT"),
      hydromad_warning = c(
        "DATA appears to cover less than 60 days.*"
      )
    ),
    `POSIXct hourly` = list(
      date = as.POSIXct(sprintf("2020-01-01 %0d:00", 1:10), tz = "GMT"), hydromad_warning = c(
        "DATA appears to cover less than 60 days.*"
      )
    ),
    `POSIXlt` = list(
      date = as.POSIXlt(sprintf("2020-01-%d 3:00", 1:10), tz = "GMT"),
      hydromad_warning = c(
        "POSIXlt index converted with as.chron",
        "DATA appears to cover less than 60 days.*"
      )
    ),
    `chron` = list(
      date = chron::as.chron(sprintf("2020-01-%d 3:00", 1:10)),
      hydromad_warning = c(
        "DATA appears to cover less than 60 days.*"
      )
    ),
    `yearmon` = list(date = as.yearmon(2000 + seq(0, 9) / 12))
  )
)
