
library(hydromad)
data(Queanbeyan)

ts90s <- window(Queanbeyan, start = "1990-01-01", end = "1999-12-31")

dbmspec <- fitDbmToPeaks(hydromad(ts90s, sma = "dbm"))

dbmords <-
  tryModelOrders(update(dbmspec, routing = "armax", rfit = "sriv"),
    n = 0:3, m = 0:2, delay = 1
  )

summary_dbmords <- summary(dbmords, stats = c("ARPE", "r.squared", "r.sq.log"))

save(summary_dbmords, file = "ihacres_tf-orders.Rdata")
