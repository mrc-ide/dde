## This needs to be run *before* running the model to avoid the solver
## unlocking problem in deSolve.
loadNamespace("deSolve")

Sys.setenv(DDE_INCLUDE = system.file("include", package = "dde"))

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

run_seir_deSolve <- function(times, tol = 1e-7) {
  b <- 1 / 10
  N <- 1e7
  beta <- 10
  sigma <- 1 / 3
  delta <- 1 / 21
  lat_hum <- 14
  I0 <- 1

  Births <- N*b
  ## i.e. proportion of humans surviving the latent period
  surv <- exp(-b*lat_hum)

  t0 <- NULL
  y0 <- NULL
  lag <- NULL

  initial <- function(t=0, pars=NULL) {
    if ("I0" %in% names(pars)) {
      I0 <<- pars$I0
    }
    t0 <<- t
    y0 <<- c(S = N - I0,E = 0, I = I0, R = 0)
    lag <<- make_lagvalue(t0, y0)
    y0
  }

  derivs <- function(t, y, .) {
    S <- y[[1L]]
    E <- y[[2L]]
    I <- y[[3L]]
    R <- y[[4L]]

    ## people developing latent infection
    new_inf <- beta*S*I/N

    ## people that become latent 'lat_hum' days ago, less those that
    ## died during that time
    S_lag <- lag(t, lat_hum, 1L)
    I_lag <- lag(t, lat_hum, 3L)
    lag_inf <- S_lag * I_lag * beta * surv / N

    dS <- Births  - b * S - new_inf + delta * R
    dE <- new_inf - lag_inf - b * E
    dI <- lag_inf - (b + sigma) * I
    dR <- sigma * I - b * R - delta * R

    list(c(dS, dE, dI, dR))
    ## Output variables
    ## c(prev = I/N, Hpop = S+E+I+R))
  }

  make_lagvalue <- function(t0, y0) {
    force(t0)
    y0 <- unname(y0)
    function(t, lag, nr=0L) {
      t1 <- t - lag
      if (t1 > t0) {
        deSolve::lagvalue(t1, nr)
      } else if (length(nr) == 1 && nr == 0L) {
        y0
      } else {
        y0[nr]
      }
    }
  }

  deSolve::dede(initial(), times, derivs,
                atol = tol, rtol = tol)[, -1, drop=FALSE]
}

run_lorenz_dde <- function(times, tol = 1e-7, n_history = 0L) {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  dopri5(y, times, "lorenz", p, atol = tol, rtol = tol, n_history = n_history,
         dllname = "lorenz")
}

run_seir_dde <- function(times, tol = 1e-7, keep_history = FALSE) {
  n_history <- 1000L
  p <- numeric(0)
  y0 <- c(1e7 - 1, 0, 1, 0)
  dopri5(y0, times, "seir", p, atol = tol, rtol = tol, n_history = n_history,
         dllname = "seir", keep_history = keep_history)
}

compile_shlib <- function(path) {
  Sys.setenv("R_TESTS" = "")
  R <- file.path(R.home(), "bin", "R")
  message("Compiling ", path)
  ok <- system2(R, c("CMD", "SHLIB", path), stdout=FALSE, stderr=FALSE)
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