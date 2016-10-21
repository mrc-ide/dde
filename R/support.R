dde_example_path <- function() {
  path <- system.file("examples", package = "dde")
  if (!nzchar(path)) {
    ## When running under at least *some* versions of
    ## devtools/testthat this is needed (but not all I think).  This
    ## is required to simplify using files from within the package in
    ## tests.  Related to workarounds in helper-dde for trying to find
    ## the include directory.  Bizarrely, I can't get this to reliably
    ## fail.
    path <- system.file("inst/examples", package = "dde") # nocov
  }
  path
}
