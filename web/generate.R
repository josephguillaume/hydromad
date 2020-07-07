## Copyright (C) Felix Andrews <felix@nfrac.org>


webitem.default <-
  function(name, ...,
           do.example = TRUE, do.helplink = TRUE,
           content = NULL) {
    c(
      webitem.basic(name, ...),
      if (do.helplink) webitem.helplink(name, ...),
      content,
      if (do.example) webitem.lattice.example(name, ...),
      webitem.codelink(name, ...)
    )
  }

webitem.basic <- function(name, ..., desc = NULL) {
  c(
    '  <h2 class="itemname">', name, "  </h2>",
    '  <div class="itemdesc">', desc, "  </div>"
  )
}

## runs Rdconv and puts HTML into 'man' dir.
## 'man.src.dir': path to man/.Rd files
webitem.helplink <-
  function(name, ..., package,
           helpname = name, examplename = helpname,
           man.src.dir, Links, Links2 = character(),
           helplinktext = "Usage, Details, Examples") {
    ## generate HTML man page file
    if (!file.exists("man")) {
      dir.create("man")
    }
    manhtml <- paste("man/", helpname, ".html", sep = "")
    manRd <- paste(man.src.dir, helpname, ".Rd", sep = "")
    tools::Rd2HTML(manRd,
      out = manhtml, package = package,
      Links = Links, Links2 = Links2
    )
    message(manhtml, " generated")
    ## generated HTML is invalid; fix it:
    tmp <- readLines(manhtml)
    # tmp <- gsub('</p>\n<p>', '<br/>', tmp)
    # tmp <- gsub('</?p>', '', tmp)
    tmp <- sub("^<!DOCTYPE .*$", "", tmp)
    tmp <- sub("^<meta .*$", "", tmp)
    tmp <- sub("^<link .*$", "", tmp)
    tmp <- gsub("<hr/?>", "", tmp)
    write(tmp, manhtml)
    ## generate HTML content
    c(
      "  <p>",
      sprintf(
        '  <a href="man/%s.html" class="helplink">',
        helpname
      ),
      helplinktext, "</a>",
      "  </p>"
    )
  }

webitem.codelink <-
  function(name, ...,
           helpname = name,
           codefile = paste(helpname, ".R", sep = ""),
           code.url = NA) {
    if (is.na(code.url) || is.na(codefile)) {
      return("")
    }
    url <- sprintf(code.url, codefile)
    sprintf(
      '<p><a href="%s" class="codelink">Source code</a></p>',
      url
    )
  }

## runs examples and puts images into 'plots' dir.
##   (currently uses dev2bitmap to make pngs)
webitem.lattice.example <-
  function(name, plotnumber = 1, ..., package,
           helpname = name, examplename = helpname,
           themes = list(default = list(theme = standard.theme("pdf"))),
           width = 500, height = 350, rerun = FALSE,
           image.src.base = "",
           call.width = 48, call.maxchar = 250) {
    ## for filenames and DOM ids
    okname <- gsub(" ", "_", name)

    if (plotnumber == 0) {
      return("")
    }

    ## we want to be able to run example() for each function
    ## but only to keep *one* of the lattice plots produced
    ## (specified by number)

    ## the following approach will only work for examples which
    ## don't include post-plotting annotations, or grid.new etc.

    ## keep track of examples as they run
    tracker <- new.env()
    tracker$plots <- list()
    tracker$counter <- 0

    OPAR <- trellis.par.get()

    ## set the lattice print function to store the target plot.
    lattice.options(print.function = function(x, ...) {
      plot(x, ...)
      tracker$counter <- tracker$counter + 1
      tracker$plots[[tracker$counter]] <- x
    })

    ## generate PNG image of target plot in example(examplename)
    firstrun <- TRUE
    plotobj <- NULL
    pdf("Rplots.pdf")
    for (themeNm in names(themes)) {
      if (!file.exists(file.path("plots", themeNm))) {
        dir.create(file.path("plots", themeNm), recursive = TRUE)
      }
      thisfile <- file.path("plots", themeNm, paste(okname, ".png", sep = ""))
      if (firstrun) {
        filename <- thisfile
      }
      stopifnot(is.list(themes[[themeNm]]$theme))
      theme <- themes[[themeNm]]$theme
      opts <- themes[[themeNm]]$options
      trellis.par.set(theme)
      OOPT <- lattice.options(opts)
      ## need to run the whole example block to regenerate plot
      ## if there are options included in theme (and on first time etc).
      ## otherwise, can skip this and just print the trellis object.
      if (firstrun || rerun || !is.null(opts)) {
        tracker$plots <- list()
        tracker$counter <- 0
        ## run the example()s for this function
        eval.parent(call("example", examplename,
          package = package,
          local = FALSE, ask = FALSE
        ))
        if (length(tracker$plots) == 0) {
          stop("no lattice plots found in example(", name, ")")
        }
        ## plotnumber can be negative, counts back from the end
        n <- if (plotnumber >= 0) {
          plotnumber
        } else {
          tracker$counter - (-plotnumber) + 1
        }
        plotobj <- tracker$plots[[n]]
        firstrun <- FALSE
      }
      ## need to use dev2bitmap(method="pdf"); bitmap() uses postscript
      ## (using postscipt loses translucency and family="serif" fails)
      pdf("tmp.pdf", width = width / 72, height = height / 72)
      dev.control(displaylist = "enable")
      trellis.par.set(theme)
      plot(plotobj)
      dev2bitmap(thisfile,
        width = width, height = height,
        units = "px", taa = 4, gaa = 4, method = "pdf"
      )
      dev.off()
      message(thisfile, " generated")
      trellis.par.set(OPAR)
      lattice.options(OOPT)
    }
    dev.off()

    ## reset to normal plotting
    lattice.options(print.function = NULL)

    fileurl <- paste(image.src.base, filename, sep = "")

    theCall <- plotobj$call
    if (identical(theCall[[1]], quote(update)) &&
      (length(theCall) == 2)) {
      ## redundant `update` wrapper; remove for clarity
      theCall <- theCall[[2]]
    }
    itemCode <-
      toString(paste(deparse(theCall, width = call.width, control = c()),
        collapse = "\n"
      ),
      width = call.maxchar
      )

    c(
      "<p>One example:</p>",
      sprintf(
        '<img src="%s" alt="%s" width="%g" height="%g"/>',
        fileurl, name, width, height
      ),
      '<pre class="itemcode">', itemCode,
      "</pre>"
    )
  }

## reads template.html from current directory.
## replacing strings @CONTENT @NAV @INDEX @VERSIONTAG
## [over]writes index.html.
## 'spec': a nested list defining groups of items for website.
## 'webitemfun' is called for each item in each group.
generateWebsite <-
  function(package, spec, ...,
           webitemfun = webitem.default,
           toplevelcontent = NULL,
           topleveljs = NULL) {
    ## get names, aliases and descriptions from help pages
    info <- readRDS(system.file("Meta", "Rd.rds", package = package))

    ## work out which help page (which element of 'info') each item belongs to
    spec <- lapply(spec, lapply, function(x) {
      name <- x[[1]]
      helpname <- if (is.null(x$helpname)) name else x$helpname
      i <- which(sapply(info$Aliases, function(aa) helpname %in% aa))
      stopifnot(length(i) <= 1)
      if (length(i) == 0) i <- NA
      x$info.i <- i
      ## fill in 'desc' element from Title field if missing
      if (is.null(x$desc)) {
        x$desc <- info$Title[i]
        if (is.null(x$desc)) stop("no description found for ", name)
      }
      x
    })
    infoIndex <- unlist(lapply(spec, lapply, function(x) x$info.i))
    itemNames <- unlist(lapply(spec, lapply, head, 1))
    ok <- !is.na(infoIndex)
    info$itemName <- NA
    info$itemName[infoIndex[ok]] <- itemNames[ok]
    ## TODO: insert website items which do not match a help page name?
    # itemNames[!ok]
    ## remove help pages for which there is no item in spec
    info <- info[!is.na(info$itemName), ]

    ## construct local HTML links
    lens <- sapply(info$Aliases, length)
    Links <- structure(paste("#", rep.int(info$itemName, lens), sep = ""),
      names = unlist(info$Aliases)
    )

    ## set up HTML buffers
    out <- character() ## @CONTENT
    nav <- character() ## @NAV

    ## process the given spec
    commonArgs <- list(..., package = package, Links = Links)
    for (i in seq_along(spec)) {
      groupName <- names(spec)[i]
      group <- spec[[i]]
      ## group level headers / containers
      out <- c(
        out,
        "",
        sprintf('<h1 class="groupname">%s</h1>', groupName),
        ""
      )
      nav <- c(
        nav,
        sprintf(
          '<li class="navhead"><a href="#">%s</a></li>',
          groupName
        ),
        '<li class="navgroup"><ul>'
      )
      ## each item:
      for (j in seq_along(group)) {
        item <- group[[j]]
        name <- item[[1]]
        okname <- gsub(" ", "_", name)
        ## generate HTML nav entries
        navid <- paste("nav_", okname, sep = "")
        nav <- c(
          nav,
          "<li>",
          sprintf(
            '<a class="navitem" href="#%s" title="%s" id="%s">%s</a>',
            okname, item$desc, navid, name
          ),
          "</li>"
        )
        ## generate HTML content (from webitemfun)
        common.masked <- names(commonArgs) %in% names(item)
        webitemArgs <- c(
          item,
          commonArgs[!common.masked]
        )
        itemContent <- do.call(webitemfun, webitemArgs)
        out <- c(
          out,
          "",
          sprintf('<div class="item" id="%s">', okname),
          itemContent,
          "</div>", ""
        )
      }
      nav <- c(nav, "</ul></li>", "")
      out <- c(out, "", "")
    }

    ## custom top-level stuff
    out <- c(out, toplevelcontent)
    if (!is.null(topleveljs)) {
      out <- c(
        out,
        '<script type="text/javascript">',
        topleveljs,
        "</script>"
      )
    }

    ## make @INDEX
    tmp <- lapply(spec, lapply, function(x) {
      c(name = x[[1]], desc = x$desc)
    })
    idxmat <- do.call("rbind", unlist(tmp, recursive = FALSE))
    idxmat <- idxmat[order(idxmat[, 1]), ]
    index <-
      with(
        as.data.frame(idxmat),
        paste("<table>",
          paste('<tr><td><a href="', paste("#", name, sep = ""), '">',
            name, "</a></td>",
            "<td>", desc, "</td></tr>",
            sep = "", collapse = "\n"
          ),
          "</table>",
          sep = "\n"
        )
      )

    ## make @VERSIONTAG
    Rvstring <- paste(
      "R version",
      paste(R.version[c("major", "minor")], collapse = ".")
    )
    nowstring <- paste("(at ", Sys.Date(), ")", sep = "")
    vTag <- paste(
      package, "version",
      packageDescription(package)$Version,
      "on", Rvstring, nowstring
    )

    ## merge content and nav into template.html
    html <- readLines("template.html")
    html <- sub("@CONTENT", paste(out, collapse = "\n"), html)
    html <- sub("@NAV", paste(nav, collapse = "\n"), html)
    html <- sub("@INDEX", index, html)
    html <- sub("@VERSIONTAG", vTag, html)

    write(html, file = "index.html")
    message("index.html generated")
  }
