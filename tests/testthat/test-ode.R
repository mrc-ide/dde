context("ode")

test_that("ode interface", {
  tt <- seq(0, 1, length.out = 200)
  m1 <- run_lorenz_deSolve(tt)
  dimnames(m1) <- NULL

  m2 <- run_lorenz_dde(tt)
  expect_equal(t(m2), m1[-1,], tolerance=1e-6)

  m3 <- run_lorenz_dde(tt, by_column = TRUE)
  expect_equal(m3, m1[-1,], tolerance=1e-6)
  expect_identical(m3, t(m2))
})

test_that("dense output", {
  tt <- seq(0, 1, length.out = 200)
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

  ## Check column output:
  m5 <- run_lorenz_dde(tt, n_history = 1000L, by_column = TRUE)
  expect_identical(attr(m5, "history"), attr(m2, "history"))
  attr(m5, "history") <- NULL
  expect_identical(m5, t(m4))
})

test_that("output", {
  tt <- seq(0, 1, length.out = 200)

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

  res3 <- dopri5(y0, tt, "lorenz", p, dllname = "lorenz",
                 n_out = 2L, output = "lorenz_output", by_column = TRUE)
  expect_equal(attr(res3, "output"), t(output))
  attr(res3, "output") <- NULL
  expect_identical(res3, t(res2))
})

test_that("keep initial", {
  tt <- seq(0, 1, length.out = 200)

  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)

  res1 <- dopri5(y0, tt, "lorenz", p, dllname = "lorenz",
                 return_initial = TRUE)
  res2 <- dopri5(y0, tt, "lorenz", p, dllname = "lorenz",
                 return_initial = FALSE)
  expect_equal(ncol(res1), length(tt))
  expect_identical(res1[, 1], y0)
  expect_identical(res1[, -1], res2)

  res3 <- dopri5(y0, tt, "lorenz", p, dllname = "lorenz",
                 n_out = 2L, output = "lorenz_output", return_initial = TRUE)
  res4 <- dopri5(y0, tt, "lorenz", p, dllname = "lorenz",
                 n_out = 2L, output = "lorenz_output", return_initial = FALSE)
  expect_equal(ncol(res3), length(tt))
  expect_identical(res3[, 1], y0)

  out3 <- attr(res3, "output")
  out4 <- attr(res4, "output")
  expect_equal(ncol(out3), length(tt))
  expect_equal(out3[,1], range(y0)) # based on definition of output function
  expect_equal(out3[, -1], out4)

  attr(res3, "output") <- attr(res4, "output") <- NULL
  expect_equal(ncol(res3), length(tt))
  expect_identical(res3[, 1], y0)
  expect_identical(res3[, -1], res4)
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
  lorenz_output <- function(t, y, p) {
    c(min(y), max(y))
  }

  tt <- seq(0, 1, length.out = 200)
  res1 <- dopri5(y0, tt, "lorenz", p, dllname = "lorenz")
  res2 <- dopri5(y0, tt, lorenz, p)
  expect_identical(res1, res2)

  res3 <- dopri5(y0, tt, "lorenz", p, dllname = "lorenz",
                 n_out = 2L, output = "lorenz_output")
  res4 <- dopri5(y0, tt, lorenz, p, n_out = 2L, output = lorenz_output)
  expect_identical(res3, res4)
})

test_that("critical times", {
  target <- function(t, y, p) {
    if (t <= 1) {
      y
    } else {
      -5 * y
    }
  }
  tt <- seq(0, 2, length.out = 200)
  res1 <- dopri5(1, tt, target, numeric(0), return_statistics = TRUE)
  res2 <- dopri5(1, tt, target, numeric(0), tcrit=1, return_statistics = TRUE)

  s1 <- attr(res1, "statistics")
  s2 <- attr(res2, "statistics")
  expect_lt(s2[["n_step"]], s1[["n_step"]])
})

test_that("names", {
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
  lorenz_output <- function(t, y, p) {
    c(min(y), max(y))
  }

  tt <- seq(0, 1, length.out = 200)

  nms <- letters[1:3]
  cmp <- list(nms, NULL)

  expect_null(dimnames(dopri5(y0, tt, lorenz, p)))
  expect_equal(dimnames(dopri5(y0, tt, lorenz, p, ynames=nms)), cmp)
  expect_equal(dimnames(dopri5(setNames(y0, nms), tt, lorenz, p)), cmp)
  expect_null(dimnames(dopri5(setNames(y0, nms), tt, lorenz, p, ynames=FALSE)))

  expect_error(dopri5(y0, tt, lorenz, p, ynames=nms[1]),
               "ynames must be the same length as y")
  expect_error(dopri5(y0, tt, lorenz, p, ynames=1),
               "Invalid value for ynames")

  ## Similar for output names:
  onms <- LETTERS[1:2]
  ocmp <- list(onms, NULL)
  expect_null(dimnames(attr(dopri5(y0, tt, lorenz, p, n_out = 2L,
                                   output = lorenz_output), "output")))
  expect_null(dimnames(attr(dopri5(y0, tt, lorenz, p, n_out = 2L,
                                   output = lorenz_output, outnames = NULL),
                            "output")))
  expect_equal(dimnames(attr(dopri5(y0, tt, lorenz, p, n_out = 2L,
                                    output = lorenz_output, outnames = onms),
                             "output")), ocmp)
  expect_error(dopri5(y0, tt, lorenz, p, n_out = 2L,
                      output = lorenz_output, outnames = nms),
               "outnames must have length n_out")
  expect_error(dopri5(y0, tt, lorenz, p, n_out = 2L,
                      output = lorenz_output, outnames = 1),
               "Invalid value for outnames")

  ## Check both together:
  res <- dopri5(y0, tt, lorenz, p, n_out = 2L,
                output = lorenz_output, ynames = nms, outnames = onms)
  expect_equal(dimnames(res), cmp)
  expect_equal(dimnames(attr(res, "output")), ocmp)
})

test_that("return_time", {
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

  tt <- seq(0, 1, length.out = 200)

  expect_equal(nrow(dopri5(y0, tt, lorenz, p)), 3)
  res <- dopri5(y0, tt, lorenz, p, return_time = TRUE)
  expect_equal(nrow(res), 4)
  expect_equal(res[1, ], tt[-1L])

  ## include the first time:
  expect_equal(
    dopri5(y0, tt, lorenz, p, return_time = TRUE, return_initial = TRUE)[1, ],
    tt)

  ## with names:
  nms <- letters[1:3]
  res <- dopri5(y0, tt, lorenz, p, return_time = TRUE, ynames = nms)
  expect_equal(rownames(res), c("time", nms))

  res <- dopri5(y0, tt, lorenz, p, return_time = TRUE, ynames = nms,
                by_column = TRUE)
  expect_equal(colnames(res), c("time", nms))
})

test_that("deSolve mode", {
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

  tt <- seq(0, 1, length.out = 200)

  res <- dopri5(y0, tt, lorenz, p, deSolve_compatible = TRUE)
  expect_equal(dim(res), c(length(tt), 4))
  expect_equal(dimnames(res), list(NULL, c("time", 1:3)))
  expect_equal(res[, 1], tt)
})
