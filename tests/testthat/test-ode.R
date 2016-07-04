context("ode")

test_that("ode interface", {
  tt <- seq(0, 1, length.out=200)
  m1 <- run_lorenz_deSolve(tt)
  dimnames(m1) <- NULL

  m2 <- run_lorenz_dde(tt)

  expect_equal(m1[-1,], m2, tolerance=1e-6)
})

test_that("dense output", {
  tt <- seq(0, 1, length.out=200)
  m1 <- run_lorenz_deSolve(tt)
  dimnames(m1) <- NULL

  m2 <- run_lorenz_dde(c(0, max(tt)), n_history = 1000L)
  ## single row of output:
  expect_equal(dim(m2), c(1, 3))
  h2 <- attr(m2, "history")
  expect_equal(nrow(h2), 17L) # 5 * 3 + 2

  m3 <- dopri5_interpolate(h2, tt)
  m4 <- run_lorenz_dde(tt)

  expect_identical(m3[-1,], m4)
  expect_equal(m3, m1, tolerance=1e-6)
})
