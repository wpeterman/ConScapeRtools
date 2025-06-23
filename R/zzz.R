.onLoad <- function(libname, pkgname) {
  registerS3method("plot", "ConScapeResults", plot.ConScapeResults)
}
