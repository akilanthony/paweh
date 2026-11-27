# Declare dependencies hidden behind paweh's namespace to the Shinylive exporter.
if (FALSE) {
  library(bslib)
  library(ggplot2)
  library(mvtnorm)
  library(plotly)
  library(rlang)
  library(shiny)
  library(withr)
}

local({
  # The local test supplies the installed, pure-R package as an app file.
  package_library <- tempfile("paweh-library-")
  dir.create(package_library)
  utils::untar("paweh-library.tar.gz", exdir = package_library, tar = "internal")
  .libPaths(c(package_library, .libPaths()))

  # The package is already bundled, so dependency discovery should only
  # fetch its dependencies. Use the public constructor as-is.
  getExportedValue("paweh", "paweh_app")()
})
