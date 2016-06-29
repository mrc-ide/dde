run_lorenz <- function(times) {
  ret <- .Call("run_lorenz", times, PACKAGE="dde")
  ret[[1L]] <- t(ret[[1L]])
  ret
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
