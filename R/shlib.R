## helper function to compile shared libraries - used in the vignettes
## and tests.
shlib <- function(filename, prefix) {
  loaded <- names(getLoadedDLLs())
  src <- basename(filename)
  base <- paste0(prefix, tools::file_path_sans_ext(src))
  if (base %in% loaded) {
    return(NULL)
  }

  path <- tempfile()
  dir.create(path, FALSE, TRUE)
  file.copy(filename, path)

  makevars <- sprintf("PKG_CFLAGS = -I%s -g",
                      system.file("include", package = "dde"))
  writeLines(makevars, file.path(path, "Makevars"))

  Sys.setenv(R_TESTS = "")
  owd <- setwd(path)
  on.exit(setwd(owd))

  dll <- paste0(base, .Platform$dynlib.ext)
  R <- file.path(R.home(), "bin", "R")
  args <- c("CMD", "SHLIB", src, "-o", dll, "--preclean", "--clean")

  output <- suppressWarnings(system2(R, args, stdout = TRUE, stderr = TRUE))

  code <- attr(output, "status")
  error <- !is.null(code) && code != 0L
  if (error) {
    cat(paste0(output, "\n", collapse = ""))
    stop("Error compiling source; see above for details")
  }

  dyn.load(dll)

  list(base = base, dll = normalizePath(dll, mustWork = TRUE), path = path)
}
