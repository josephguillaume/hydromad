## generate the web site

## setwd("~/devel/hydromad/web")
## Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.06/bin/gswin64c.exe")
## source("genHydromad.R", echo = TRUE)

library(hydromad)

spec <- list()

spec[["modelling framework"]] <-
    list(
         list("hydromad"),
         list("methods...", helpname = "hydromad.object",
              codefile = "update.hydromad.R"),
         list("predict", helpname = "predict.hydromad"),
         list("simulate", helpname = "simulate.hydromad", do.example = TRUE),
         list("runlist"),
         list("hydromad.options", codefile = "options.R")
         )

spec[["assessment"]] <-
    list(
         list("summary", helpname = "summary.hydromad"),
         list("hydromad.stats"),
         list("nseStat"),
         list("nseVarTd"),
         list("xyplot", helpname = "xyplot.hydromad", do.example = TRUE),
         list("xyplot.runlist", do.example = TRUE),
         list("qqmath", -2, helpname = "xyplot.hydromad", do.example = TRUE)
         )

spec[["calibration"]] <-
    list(
         list("buildTsObjective"),
         list("objFunVal"),
         list("fitBySampling"),
         list("fitByOptim"),
         list("fitBySCE"),
         list("fitByDE"),
         list("fitByDream"),
         list("fitByCMAES"),
         list("fitByDDS"),
         list("fitByNsga2"),
         list("optimtrace", -1, do.example = TRUE, examplename = "fitByOptim"),
         list("defineFeasibleSet", do.example = TRUE)
         )

spec[["multimodel"]] <-
  list(
    list("hydromad_parallelisation",do.example=FALSE),
    list("paretoObjectivesNsga2",do.example=FALSE),
    list("paretoObjectivesVaryWeights",do.example=FALSE),
    list("crossValidate",do.example=FALSE),
    list("splitData",do.example=FALSE),
    list("paretoTimeAnalysis",do.example=FALSE)
  )

spec[["identifiability"]] <-
  list(
      list("hydromad_sensitivity",do.example=FALSE),
      list("tellTS",do.example=FALSE),
      list("runRSM",do.example=FALSE),
      list("eigen.plot",do.example=FALSE),
      list("findUnivariateBounds",do.example=FALSE)
    )
  
spec[["discrete events"]] <-
    list(
         list("eventseq", do.example = TRUE),
         list("event.xyplot", do.example = TRUE),
         list("event.xyplot.hydromad", do.example = TRUE)
         )

spec[["soil moisture accounting"]] <-
    list(
         list("cmd", helpname = "IHACRES.CMD.model", codefile = "cmd.R"),
         list("cwi", helpname = "IHACRES.CWI.model", codefile = "cwi.R"),
         list("gr4j"),
         list("awbm"),
         list("bucket"),
         list("sacramento"),
         list("snow"),
         list("scalar"),
         list("intensity"),
         list("runoffratio"),
         list("dbm")
         )

spec[["routing"]] <-
    list(
         list("armax", codefile = "armax.sim.R"),
         list("expuh", codefile = "expuh.sim.R"),
         list("lambda", codefile = "lambda.sim.R"),
         list("powuh", codefile = "powuh.sim.R"),
         list("leakyExpStore",codefile="leakyExpStore.R"),
         list("expuh3s",codefile="expuh3s.sim.R")
         ##list("maexpuh", codefile = "maexpuh.R")
         )

spec[["routing fitting"]] <-
    list(
         list("armax.ls.fit"),
         list("armax.sriv.fit"),
         list("armax.inverse.fit"),
         #"lambda.inverse.fit"
         list("armax.inverse.sim", do.example = TRUE),
         list("tryModelOrders"),
         list("estimateDelay"),
         list("estimateDelayFrac"),
         list("deconvolution.uh")
         )

spec[["utilities"]] <-
    list(
         list("convertFlow"),
         list("rollccf", do.example = TRUE),
         list("parameterSets"),
         list("evalPars"),
         ##list("getFreeParsRanges"), ##same page as evalPars
         list("evalParsRollapply"),
         list("rotatedSampling"),
         list("gr4j.transformpar"),
         list("observed"),
         list("SCEoptim", codefile = "sce.R")
         )

spec[["datasets"]] <-
    list(
         list("BinghamTrib"),
         list("Canning"),
         list("Cotter"),
         #list("Molonglo"),
         list("Murrindindi"),
         list("Queanbeyan"),
         list("SalmonBrook"),
         list("Wye"),
         list("HydroTestData")
         )

spec[["datasets"]] <- lapply(spec[["datasets"]], function(x) {
    x$do.example <- TRUE
    x$codefile <- paste("../data/", x[[1]], ".R", sep = "")
    x
})

#spec[["wetlands"]] <-
#    list(
#         list("swimp"),
#         list("poweroid")
#         )

#source("http://latticeextra.r-forge.r-project.org/generate.R")
source("generate.R")

## stop on errors
lattice.options(panel.error = NULL)

ltheme <- custom.theme.2()
ltheme$strip.background$col <- grey(7/8)
ltheme$strip.shingle$col <- grey(6/8)

imageSrcBase <- ""

generateWebsite("hydromad", spec = spec, do.example = FALSE, ## default
                man.src.dir = "../man/",
                image.src.base = imageSrcBase,
                topleveljs = paste('var imageSrcBase = "', imageSrcBase, '";', sep = ""),
                code.url = "http://github.com/josephguillaume/hydromad/tree/master/R/%s",
                themes = list(custom_theme_2 = list(theme = ltheme)))
