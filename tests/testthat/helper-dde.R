## This needs to be run *before* running the model to avoid the solver
## unlocking problem in deSolve.
loadNamespace("deSolve")

## A simple Lorenz attractor
run_lorenz_deSolve <- function(times, tol = 1e-7) {
  sigma <- 10.0
  R     <- 28.0
  b     <-  8.0 / 3.0

  initial <- function(t=0, pars=NULL) {
    c(10, 1, 1)
  }

  derivs <- function(t, y, .) {
    y1 <- y[[1L]]
    y2 <- y[[2L]]
    y3 <- y[[3L]]
    list(c(sigma * (y2 - y1),
           R * y1 - y2 - y1 * y3,
           -b * y3 + y1 * y2))
  }

  deSolve::ode(initial(), times, derivs,
               atol = tol, rtol = tol)[, -1, drop=FALSE]
}

run_lorenz_dde <- function(times, tol = 1e-7, n_history = 0L) {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  dopri5(y, times, "lorenz", p, atol = tol, rtol = tol, n_history = n_history,
         dllname = "lorenz")
}

compile_shlib <- function(path) {
  Sys.setenv("R_TESTS" = "")
  R <- file.path(R.home(), "bin", "R")
  message("Compiling ", path)
  ok <- system2(R, c("CMD", "SHLIB", path), stdout=FALSE, stderr=FALSE)
  shlib <- sub("\\.c$", .Platform$dynlib.ext, path)
  base <- sub("\\.c$", "", basename(path))
  if (base %in% names(getLoadedDLLs())) {
    dyn.unload(shlib)
  }
  dyn.load(shlib)
}

cleanup_objects <- function() {
  files <- dir(pattern="\\.(o|so|dll)$")
  if (length(files) > 0L) {
    file.remove(files)
  }
  invisible()
}

prepare_all <- function() {
  cleanup_objects()
  files <- dir(pattern = "\\.c$")
  for (f in files) {
    compile_shlib(f)
  }
}

prepare_all()
