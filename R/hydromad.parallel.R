hydromad.parallel <- function(settings) {
  output <- list(
    method = "none",
    export = c(),
    packages = "hydromad",
    async = FALSE
  )
  if (missing(settings)) {
    return(output)
  }
  ## Backwards compatible-case, specifying only name of method
  if (is.character(settings)) {
    return(modifyList(output, list(method = settings)))
  }
  if (is.list(settings)) {
    return(modifyList(output, settings))
  }
  return(output)
}
