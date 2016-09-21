## I'm going to need this in both the tests and within the vignette,
## so this needs to be defined here.  I'll not expose it for now
## though.  It would be nice if this was something available via
## package:tools or something similar though.
compile_shlib <- function(path, verbose = interactive()) {
  stdout <- if (verbose) '' else FALSE
  Sys.setenv("R_TESTS" = "")
  R <- file.path(R.home(), "bin", "R")
  owd <- setwd(dirname(path))
  on.exit(setwd(owd))

  if (verbose) {
    message("Compiling ", path)
  }
  ok <- system2(R, c("CMD", "SHLIB", path), stdout=stdout, stderr=stdout)
  if (ok != 0L) {
    stop(sprintf("compilation of %s failed", path))
  }
  shlib <- sub("\\.c$", .Platform$dynlib.ext, basename(path))
  base <- sub("\\.c$", "", basename(path))
  if (base %in% names(getLoadedDLLs())) {
    dyn.unload(shlib)
  }
  dyn.load(shlib)
}

dde_example_path <- function() {
  path <- system.file("examples", package = "dde")
  if (!nzchar(path)) {
    ## When running under at least *some* versions of
    ## devtools/testthat this is needed (but not all I think).  This
    ## is required to simplify using files from within the package in
    ## tests.  Related to workarounds in helper-dde for trying to find
    ## the include directory.  Bizarrely, I can't get this to reliably
    ## fail.
    path <- system.file("inst/examples", package = "dde")
  }
  path
}
