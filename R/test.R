run_lorenz_dde <- function(times, n_history = 100L) {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  dopri5(y, tt, "lorenz", p, atol=1e-7, rtol=1e-7, n_history=n_history)
}

run_lorenz_deSolve <- function(times) {
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

  deSolve::ode(initial(), times, derivs, atol=1e-7, rtol=1e-7)[, -1, drop=FALSE]
}
