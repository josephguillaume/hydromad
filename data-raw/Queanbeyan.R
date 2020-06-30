Queanbeyan <- local({

  ## read files from doc directory
  pqdat <- read.table("Queanbeyan.dat",
    sep = ",", col.names = c("P", "Q", "Date"), as.is = TRUE
  )
  tdat <- read.table("t_cotter.csv",
    sep = ",", col.names = c("T", "Date"), as.is = TRUE
  )
  ## convert dates
  pqdat$Date <- as.Date(pqdat$Date, "%d/%m/%Y")
  tdat$Date <- as.Date(tdat$Date, "%d/%m/%Y")

  ## convert missing values
  pqdat$P[pqdat$P < 0] <- NA
  pqdat$Q[pqdat$Q < 0] <- NA
  tdat <- subset(tdat, !is.na(Date))

  ## convert from ML to mm (i.e. divide by catchment area km^2)
  pqdat$Q <- pqdat$Q / 490
  #    pqdat$Q <- convertFlow(pqdat$Q, from = "ML", area.km2 = 490)

  ## zoo objects
  library(zoo)
  tsPQ <- zoo(pqdat[, 1:2], order.by = pqdat$Date, frequency = 1)
  tsT <- zoo(tdat[, 1], order.by = tdat$Date, frequency = 1)
  tsPQE <- merge(tsPQ, E = tsT, all = FALSE)
  na.trim(tsPQE)
})
usethis::use_data(Queanbeyan, overwrite = TRUE)
