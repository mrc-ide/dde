context("ode")

test_that("ode interface", {
  tt <- seq(0, 1, length.out=200)
  m1 <- run_lorenz_deSolve(tt)
  dimnames(m1) <- NULL

  m2 <- run_lorenz_dde(tt)

  expect_equal(m1[-1,], t(m2), tolerance=1e-6)
})

test_that("dense output", {
  tt <- seq(0, 1, length.out=200)
  m1 <- run_lorenz_deSolve(tt)
  dimnames(m1) <- NULL

  m2 <- run_lorenz_dde(c(0, max(tt)), n_history = 1000L)
  ## single row of output:
  expect_equal(dim(m2), c(3, 1))
  h2 <- attr(m2, "history")
  expect_equal(nrow(h2), 17L) # 5 * 3 + 2

  m3 <- dopri5_interpolate(h2, tt)
  m4 <- run_lorenz_dde(tt)

  expect_identical(m3[-1,], t(m4))
  expect_equal(m3, m1, tolerance=1e-6)
})

test_that("output", {
  tt <- seq(0, 1, length.out=200)

  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)

  res1 <- dopri5(y0, tt, "lorenz", p, dllname = "lorenz")
  res2 <- dopri5(y0, tt, "lorenz", p, dllname = "lorenz",
                 n_out = 2L, output = "lorenz_output")

  expect_equal(names(attributes(res1)), "dim")
  output <- attr(res2, "output")
  expect_equal(dim(output), c(2L, ncol(res1)))
  expect_equal(output[1L, ], pmin(res1[1,], res1[2,], res1[3,]))
  expect_equal(output[2L, ], pmax(res1[1,], res1[2,], res1[3,]))

  attr(res2, "output") <- NULL
  expect_identical(res1, res2)
})

test_that("R interface", {
  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)

  lorenz <- function(t, y, p) {
    sigma <- p[[1L]]
    R <- p[[2L]]
    b <- p[[3L]]
    c(sigma * (y[[2L]] - y[[1L]]),
      R * y[[1L]] - y[[2L]] - y[[1L]] * y[[3L]],
      -b * y[[3L]] + y[[1L]] * y[[2L]])
  }

  tt <- seq(0, 1, length.out=200)
  res1 <- dopri5(y0, tt, "lorenz", p, dllname = "lorenz")
  res2 <- dopri5(y0, tt, lorenz, p)
  expect_identical(res1, res2)

  res2 <- dopri5(y0, tt, lorenz, p, output = lorenz_output)
})
