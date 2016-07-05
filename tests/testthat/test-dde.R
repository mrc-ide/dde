context("dde")

test_that("dde", {
  ## The pre-lag part of the integration
  tt1 <- seq(0, 14, length.out=101)
  ## Post lag, but the first bit
  tt2 <- seq(0, 30, length.out=301)
  ## The entire interesting phase
  tt3 <- seq(0, 200, length.out=301)

  yy1 <- run_seir_deSolve(tt1)
  yy2 <- run_seir_deSolve(tt2)
  yy3 <- run_seir_deSolve(tt3)
  dimnames(yy1) <- dimnames(yy2) <- dimnames(yy3) <- NULL

  ## Before the lag; this one is easy:
  zz1 <- run_seir_dde(tt1)
  expect_equal(t(zz1), yy1[-1, ])

  zz2 <- run_seir_dde(tt2)
  expect_equal(t(zz2), yy2[-1, ])

  zz3 <- run_seir_dde(tt3)
  expect_equal(t(zz3), yy3[-1, ], tolerance = 1e-6)
})

test_that("output", {
  tt <- seq(0, 200, length.out=301)
  p <- numeric(0)
  y0 <- c(1e7 - 1, 0, 1, 0)
  res1 <- dopri5(y0, tt, "seir", p, n_history = 1000L,
                 dllname = "seir", keep_history = FALSE)
  expect_equal(names(attributes(res1)), "dim")

  res2 <- dopri5(y0, tt, "seir", p, n_history = 1000L,
                 n_out = 1L, output = "seir_output",
                 dllname = "seir", keep_history = FALSE)

  output <- attr(res2, "output")
  expect_equal(dim(output), c(1L, ncol(res1)))
  expect_equal(drop(output), colSums(res2), tolerance = 1e-14)

  attr(res2, "output") <- NULL
  expect_identical(res1, res2)
})

test_that("R interface", {
  tt <- seq(0, 200, length.out=301)
  p <- numeric(0)
  y0 <- c(1e7 - 1, 0, 1, 0)
  res1 <- dopri5(y0, tt, "seir", p, n_history = 1000L,
                 dllname = "seir", keep_history = FALSE)

  seir <- function(t, y, p) {
    b <- 1.0 / 10.0
    N <- 1e7
    beta <- 10.0
    sigma <- 1.0 / 3.0
    delta <- 1.0 / 21.0
    lat_hum <- 14.0
    Births <- N * b
    surv <- exp(-b * lat_hum)

    S <- y[[1L]]
    E <- y[[2L]]
    I <- y[[3L]]
    R <- y[[4L]]

    tau <- t - lat_hum
    if (p == "one") {
      S_lag <- ylag(tau, 1L)
      I_lag <- ylag(tau, 3L)
    } else if (p == "idx") {
      y_lag <- ylag(tau, c(1L, 3L))
      S_lag <- y_lag[[1L]]
      I_lag <- y_lag[[2L]]
    } else {
      y_lag <- ylag(tau)
      S_lag <- y_lag[[1L]]
      I_lag <- y_lag[[3L]]
    }

    new_inf <- beta * S * I / N
    lag_inf <- beta * S_lag * I_lag * surv / N

    c(Births  - b * S - new_inf + delta * R,
      new_inf - lag_inf - b * E,
      lag_inf - (b + sigma) * I,
      sigma * I - b * R - delta * R)
  }

  res2 <- dopri5(y0, tt, seir, "one", n_history = 1000L, keep_history = FALSE)
  res3 <- dopri5(y0, tt, seir, "idx", n_history = 1000L, keep_history = FALSE)
  res4 <- dopri5(y0, tt, seir, "all", n_history = 1000L, keep_history = FALSE)

  expect_equal(res2, res1, tolerance = 1e-14)
  expect_identical(res3, res2)
  expect_identical(res4, res2)
})
