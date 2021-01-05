context("ode")

test_that("ode interface", {
  skip_if_not_installed("deSolve")
  tt <- seq(0, 1, length.out = 200)
  m1 <- run_lorenz_deSolve(tt)
  dimnames(m1) <- NULL

  m2 <- run_lorenz_dde(tt)
  expect_equal(m2, m1, tolerance = 1e-6)

  m3 <- run_lorenz_dde(tt, return_minimal = TRUE)
  expect_equal(t(m3), m1[-1, -1, drop = FALSE], tolerance = 1e-6)
  expect_identical(m3, t(m2[-1, -1]))
})

test_that("integer time input", {
  tt <- seq(0, 10, by = 1)
  expect_identical(run_lorenz_dde(as.integer(tt)),
                   run_lorenz_dde(tt))
})

test_that("zero time difference", {
  tt <- c(0, 0, 1, 1, 2, 2, 2, 3, 3, 3)
  res1 <- run_lorenz_dde(tt, return_initial = TRUE)
  res2 <- run_lorenz_dde(unique(tt), return_initial = TRUE)
  expect_identical(res1, res2[tapply(tt, tt), ])
})

test_that("ode, 873 stepper", {
  tt <- seq(0, 1, length.out = 200)
  m5 <- run_lorenz_dde(tt)
  m8 <- run_lorenz_dde(tt, method = "dopri853")
  expect_equal(m5, m8, tolerance = 1e-7)
  expect_false(identical(m5, m8))
})

test_that("dense output", {
  skip_if_not_installed("deSolve")
  tt <- seq(0, 1, length.out = 200)
  m1 <- run_lorenz_deSolve(tt)
  dimnames(m1) <- NULL

  for (method in dopri_methods()) {
    m2 <- run_lorenz_dde(c(0, max(tt)), n_history = 1000L, method = method)
    ## single row of output:
    expect_equal(dim(m2), c(2, 4))
    h2 <- attr(m2, "history")
    expect_is(h2, "dopri_history")
    expect_equal(nrow(h2),
                 if (method == "dopri5") 17L else 26L) # c(5, 8) * 3 + 2
    expect_identical(attr(h2, "n"), 3L)

    expect_output(print(h2), "<dopri_history>", fixed = TRUE)

    m3 <- dopri_interpolate(h2, tt)
    m4 <- run_lorenz_dde(tt, method = method)

    expect_equal(m3, m4[, -1], tolerance = 1e-10)
    expect_equal(m3, m1[, -1], tolerance = 1e-6)

    expect_equal(dopri_interpolate(m2, tt), m3)

    ## Check row output:
    m5 <- run_lorenz_dde(tt, n_history = 1000L, return_by_column = FALSE,
                         method = method)
    expect_identical(attr(m5, "history"), attr(m2, "history"))
    attr(m5, "history") <- NULL
    expect_identical(m5, t(m4))

    expect_error(dopri_interpolate(h2, tt - 1),
                 "Time falls outside of range of known history")
    expect_error(dopri_interpolate(h2, tt + 1),
                 "Time falls outside of range of known history")

    expect_error(dopri_interpolate(h2[-1, ], tt),
                 "Corrupt history object")
    h2 <- h2[-1, ]
    attr(h2, "n") <- 3L
    class(h2) <- "dopri_history"
    expect_error(dopri_interpolate(h2, tt),
                 "Corrupt history object: incorrect number of rows")
  }
})

test_that("output", {
  tt <- seq(0, 1, length.out = 200)

  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)

  res1 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz")
  res2 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz",
                return_output_with_y = FALSE,
                n_out = 2L, output = "lorenz_output")

  expect_equal(names(attributes(res1)), "dim")
  output <- attr(res2, "output")
  expect_equal(dim(output), c(nrow(res1), 2L))
  expect_equal(output[, 1L], pmin(res1[, 2], res1[, 3], res1[, 4]))
  expect_equal(output[, 2L], pmax(res1[, 2], res1[, 3], res1[, 4]))

  attr(res2, "output") <- NULL
  expect_identical(res1, res2)

  res3 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz",
                n_out = 2L, output = "lorenz_output",
                return_output_with_y = FALSE, return_by_column = FALSE)
  expect_equal(attr(res3, "output"), t(output))
  attr(res3, "output") <- NULL
  expect_identical(res3, t(res2))

  res4 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz",
                n_out = 2L, output = "lorenz_output",
                return_by_column = FALSE, return_output_with_y = TRUE)
  expect_null(attr(res4, "output"))
  expect_equal(dim(res4), c(6, length(tt)))
  expect_equal(res4[1:4, ], res3, check.attributes = FALSE)
  expect_equal(res4[5:6, ], t(output), check.attributes = FALSE)

  res5 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz",
                n_out = 2L, output = "lorenz_output",
                return_by_column = FALSE, return_output_with_y = TRUE,
                return_time = FALSE, return_initial = FALSE)
  expect_null(attr(res5, "output"))
  expect_equal(dim(res5), c(5, length(tt) - 1L))
  expect_equal(res5[1:3, ], res3[-1, -1], check.attributes = FALSE)
  expect_equal(res5[4:5, ], t(output[-1, ]), check.attributes = FALSE)
})

test_that("keep initial", {
  tt <- seq(0, 1, length.out = 200)

  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)

  res1 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz",
                return_initial = TRUE)
  res2 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz",
                return_initial = FALSE)
  expect_equal(nrow(res1), length(tt))
  expect_identical(res1[1, -1], y0)
  expect_identical(res1[-1, ], res2)

  res3 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz",
                n_out = 2L, output = "lorenz_output", return_initial = TRUE,
                return_output_with_y = FALSE)
  res4 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz",
                n_out = 2L, output = "lorenz_output", return_initial = FALSE,
                return_output_with_y = FALSE)
  expect_equal(nrow(res3), length(tt))
  expect_identical(res3[1, -1], y0)

  out3 <- attr(res3, "output")
  out4 <- attr(res4, "output")
  expect_equal(nrow(out3), length(tt))
  expect_equal(out3[1, ], range(y0)) # based on definition of output function
  expect_equal(out3[-1, ], out4)

  output <- attr(res3, "output")
  attr(res3, "output") <- attr(res4, "output") <- NULL
  expect_equal(nrow(res3), length(tt))
  expect_identical(res3[1, -1], y0)
  expect_identical(res3[-1, ], res4)

  res5 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz",
                n_out = 2L, output = "lorenz_output", return_initial = TRUE,
                return_output_with_y = TRUE)
  res6 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz",
                n_out = 2L, output = "lorenz_output", return_initial = FALSE,
                return_output_with_y = TRUE)
  expect_null(attr(res5, "output"))
  expect_null(attr(res6, "output"))
  expect_equal(res5[, 1:4], res3)
  expect_equal(res6[, 1:4], res4)
  expect_equal(res5[, 5:6], output)
  expect_equal(res6[, 5:6], output[-1, ])
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
  res1 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz")
  res2 <- dopri(y0, tt, lorenz, p)
  expect_equal(res1, res2, tolerance = 1e-10)

  res3 <- dopri(y0, tt, "lorenz", p, dllname = "dde_lorenz",
                n_out = 2L, output = "lorenz_output")
  res4 <- dopri(y0, tt, lorenz, p, n_out = 2L, output = lorenz_output)
  expect_equal(res3, res4, tolerance = 1e-10)
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
  res1 <- dopri(1, tt, target, numeric(0), return_statistics = TRUE)
  res2 <- dopri(1, tt, target, numeric(0), tcrit = 1, return_statistics = TRUE)
  res3 <- dopri(1, tt, target, numeric(0), tcrit = c(-1, 1),
                return_statistics = TRUE)

  s1 <- attr(res1, "statistics")
  s2 <- attr(res2, "statistics")
  s3 <- attr(res3, "statistics")
  expect_lt(s2[["n_step"]], s1[["n_step"]])
  expect_equal(s2, s3)
  expect_identical(res2, res3)

  expect_is(attr(res1, "step_size"), "numeric")

  ## I also need to check a pathology here:
  res4 <- dopri(1, tt, target, numeric(0), tcrit = c(tt[[1]], 1),
                return_statistics = TRUE)
  expect_equal(res4, res3)

  ## Bunch of pathalogical times:
  res5 <- dopri(1, tt, target, numeric(0), tcrit = rep(tt[[1]], 3),
                return_statistics = TRUE)
  expect_equal(res5, res1)

  ## Pile up some times in the middle:
  res6 <- dopri(1, tt, target, numeric(0), tcrit = rep(1, 3),
                return_statistics = TRUE)
  expect_equal(res4, res3)
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
  cmp <- list(NULL, c("time", nms))

  expect_null(dimnames(dopri(y0, tt, lorenz, p)))
  expect_equal(dimnames(dopri(y0, tt, lorenz, p, ynames = nms)), cmp)
  expect_equal(dimnames(dopri(setNames(y0, nms), tt, lorenz, p)), cmp)
  expect_null(dimnames(dopri(setNames(y0, nms), tt, lorenz, p, ynames = FALSE)))

  expect_error(dopri(y0, tt, lorenz, p, ynames = nms[1]),
               "ynames must be the same length as y")
  expect_error(dopri(y0, tt, lorenz, p, ynames = 1),
               "Invalid value for ynames")

  ## Similar for output names:
  onms <- LETTERS[1:2]
  ocmp <- list(NULL, onms)

  f <- function(..., return_output_with_y = FALSE) {
    dopri(y0, tt, lorenz, p, n_out = 2L, output = lorenz_output,
          return_output_with_y = return_output_with_y, ...)
  }
  expect_null(dimnames(attr(f(), "output")))
  expect_null(dimnames(attr(f(outnames = NULL), "output")))
  expect_equal(dimnames(attr(f(outnames = onms), "output")), ocmp)
  expect_error(f(outnames = nms), "outnames must have length n_out")
  expect_error(f(outnames = 1), "Invalid value for outnames")

  ## Check both together:
  res <- f(ynames = nms, outnames = onms)
  expect_equal(dimnames(res), cmp)
  expect_equal(dimnames(attr(res, "output")), ocmp)

  ## And check when output is combined
  m <- f(ynames = nms, outnames = NULL, return_output_with_y = TRUE)
  expect_equal(colnames(m), c("time", nms, rep("", length(onms))))
  m <- f(ynames = NULL, outnames = onms, return_output_with_y = TRUE)
  expect_equal(colnames(m), c("time", rep("", length(nms)), onms))
  m <- f(ynames = nms, outnames = onms, return_output_with_y = TRUE)
  expect_equal(colnames(m), c("time", nms, onms))

  ## And check time
  expect_null(dimnames(f(return_time = TRUE)))
  expect_equal(colnames(f(return_time = TRUE, ynames = nms)),
               c("time", nms))
  expect_equal(colnames(f(return_time = TRUE, return_output_with_y = TRUE,
                          ynames = nms)),
               c("time", nms, rep("", length(onms))))
  expect_equal(colnames(f(return_time = TRUE, return_output_with_y = TRUE,
                          ynames = nms, outnames = onms)),
               c("time", nms, onms))
  expect_equal(colnames(f(return_time = TRUE, return_output_with_y = TRUE,
                          ynames = NULL, outnames = onms)),
               c("time", rep("", length(nms)), onms))
  expect_null(colnames(f(return_time = TRUE, return_output_with_y = TRUE,
                         ynames = NULL, outnames = NULL)))
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

  res1 <- dopri(y0, tt, lorenz, p)
  res2 <- dopri(y0, tt, lorenz, p, return_time = FALSE)
  expect_equal(ncol(res1), 4)
  expect_equal(ncol(res2), 3)
  expect_equal(res1[, 1L], tt)
  expect_equal(res1[, -1L], res2)

  res3 <- dopri(y0, tt, lorenz, p, return_initial = FALSE)
  ## include the first time:
  expect_equal(res3[, 1], tt[-1])

  ## with names:
  nms <- letters[1:3]
  res <- dopri(y0, tt, lorenz, p, ynames = nms)
  expect_equal(colnames(res), c("time", nms))

  res <- dopri(y0, tt, lorenz, p, ynames = nms, return_by_column = FALSE)
  expect_equal(rownames(res), c("time", nms))
})

test_that("minimal mode", {
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
  names(y0) <- letters[1:3]

  res1 <- dopri(y0, tt, lorenz, p)
  res2 <- dopri(y0, tt, lorenz, p, return_minimal = TRUE)

  expect_equal(dim(res1), c(length(tt), 4))
  expect_equal(dim(res2), c(3, length(tt) - 1L))

  expect_equal(dimnames(res1), list(NULL, c("time", names(y0))))
  expect_null(dimnames(res2))

  expect_equal(res1[, 1], tt)
  expect_equal(t(unname(res1[-1, -1])), res2)
})

test_that("step tuning", {
  tt <- seq(0, 1, length.out = 200)
  T_IDX <- 16L
  m0 <- run_lorenz_dde(tt, n_history = 500, return_statistics = TRUE)
  t0 <- attr(m0, "history")[T_IDX, ]
  s0 <- attr(m0, "statistics")

  hmax <- max(diff(t0)) / 5
  m1 <- run_lorenz_dde(tt, step_size_max = hmax, n_history = 500)
  t1 <- attr(m1, "history")[T_IDX, ]
  expect_gt(length(t1), length(t0))
  expect_equal(max(diff(t1)), hmax)

  h0 <- diff(t0[1:2])
  m2 <- run_lorenz_dde(tt, n_history = 500, return_statistics = TRUE,
                       step_size_initial = h0)
  t2 <- attr(m0, "history")[T_IDX, ]
  expect_equal(t0, t2)
  expect_equal(attr(m2, "statistics")[["n_eval"]],
               s0[["n_eval"]] - 1L)

  expect_error(run_lorenz_dde(tt, step_max_n = length(t0) - 20),
               "too many steps")
  expect_error(run_lorenz_dde(tt, step_size_min = min(diff(t0)) * 1.1),
               "step size too small")
})


test_that("allow step size to get small", {
  tt <- seq(0, 1, length.out = 200)
  T_IDX <- 16L
  m0 <- run_lorenz_dde(tt, n_history = 500, return_statistics = TRUE)
  t0 <- attr(m0, "history")[T_IDX, ]
  s0 <- attr(m0, "statistics")

  step_size_min <- 0.01
  m1 <- run_lorenz_dde(tt,
                       step_size_min = step_size_min,
                       step_size_min_allow = TRUE,
                       n_history = 500,
                       return_statistics = TRUE)
  t1 <- attr(m1, "history")[T_IDX, ]
  ## This is really annoying, as we might be off by a _very_ small
  ## amount (~1e-18) due to general floating point horrors.
  expect_lt(min(diff(t1) - step_size_min), 1e-16)
  expect_equal(min(diff(t1)), step_size_min)
  s1 <- attr(m1, "statistics")
  expect_lt(s1[["n_eval"]], s0[["n_eval"]])
})


test_that("integrate function with no absolute error", {
  deriv <- function(t, y, p) {
    1
  }

  expect_equal(drop(dopri(0, c(0, 1), deriv, 0, return_minimal = TRUE)), 1)
})

test_that("Native Symbol interface", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 10, length.out = 101)
  func <- getNativeSymbolInfo("lorenz", PACKAGE = "dde_lorenz")
  y <- dopri(y, times, func, p)
  expect_equal(y, run_lorenz_dde(times, tol = 1e-6))
})

test_that("NULL pointer safety", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 10, length.out = 101)

  func <- getNativeSymbolInfo("lorenz", PACKAGE = "dde_lorenz")
  func <- make_null_pointer(func)

  expect_error(dopri(y, times, func, p), "null pointer")
  expect_error(dopri(y, times, "lorenz", p, output = func, n_out = 2L),
               "null pointer")
})

test_that("always return history when asked", {
  has_history <- function(x) {
    "history" %in% names(attributes(x))
  }
  has_statistics <- function(x) {
    "statistics" %in% names(attributes(x))
  }

  growth <- function(t, y, p) {
    y * p
  }
  y0 <- runif(1)
  r <- runif(1)
  tt <- seq(0, 2, length.out = 101)

  expect_false(has_history(dopri(y0, tt, growth, r)))
  expect_true(has_history(dopri(y0, tt, growth, r, n_history = 100)))

  expect_true(has_history(
    dopri(y0, tt, growth, r, n_history = 100, return_minimal = TRUE)))
  expect_true(has_history(
    dopri(y0, tt, growth, r, n_history = 100, return_by_column = TRUE)))
  expect_true(has_history(
    dopri(y0, tt, growth, r, n_history = 100, return_initial = TRUE)))
  expect_true(has_history(
    dopri(y0, tt, growth, r, n_history = 100, return_time = TRUE)))

  expect_false(has_statistics(dopri(y0, tt, growth, r)))
  expect_true(has_statistics(dopri(y0, tt, growth, r,
                                   return_statistics = TRUE)))

  expect_true(has_statistics(
    dopri(y0, tt, growth, r, return_statistics = TRUE,
          return_minimal = TRUE)))
  expect_true(has_statistics(
    dopri(y0, tt, growth, r, return_statistics = TRUE,
          return_by_column = TRUE)))
  expect_true(has_statistics(
    dopri(y0, tt, growth, r, return_statistics = TRUE,
          return_initial = TRUE)))
  expect_true(has_statistics(
    dopri(y0, tt, growth, r, return_statistics = TRUE,
          return_time = TRUE)))
})

test_that("negative time", {
  skip_if_not_installed("deSolve")
  growth <- function(t, y, p) {
    y * p
  }

  set.seed(1)
  y0 <- runif(4)
  r <- runif(4)

  ## True solution is:
  ##
  ##   exp(r * t) * y0
  ##
  ## which we can compute easily enough foe this problem:
  tt <- seq(0, -2, length.out = 11)
  real <- t(y0 * outer(r, tt, function(t, r) exp(r * t)))
  real_fwd <- t(y0 * outer(-r, -tt, function(t, r) exp(r * t)))
  expect_identical(real, real_fwd)

  ## Then solve the system using deSolve (sanity check for the above)
  cmp <- deSolve::ode(y0, tt, function(...) list(growth(...)), r)
  cmp_fwd <- deSolve::ode(y0, -tt, function(...) list(growth(...)), -r)
  expect_equal(unname(cmp[, -1]), real, 5e-6)
  expect_equal(cmp[, -1], cmp_fwd[, -1], tolerance = 1e-16)

  ## and dde:
  res <- dopri(y0, tt, growth, r, n_history = 100)
  res_fwd <- dopri(y0, -tt, growth, -r, n_history = 100)

  expect_equal(res[, -1], res_fwd[, -1], tolerance = 1e-16)

  ## Interestingly, we do stop at the same points here, so the error
  ## calculation and underlying stepping is probably OK.
  t_idx <- length(y0) * 5 + 1
  h1 <- attr(res, "history")
  h2 <- attr(res_fwd, "history")
  expect_equal(-h1[t_idx, ], h2[t_idx, ], tolerance = 1e-16)
  expect_equal(-h1[t_idx + 1, ], h2[t_idx + 1, ], tolerance = 1e-16)
  expect_equal(h1[1:20, ], h2[1:20, ], tolerance = 1e-16)
})

test_that("negative time with tcrit", {
  skip_if_not_installed("deSolve")
  target1 <- function(t, y, p) {
    if (t <= 1) y else -5 * y
  }
  target2 <- function(t, y, p) {
    if (t >= -1) -y else 5 * y
  }

  ## Here's our system.  It's very silly again.
  tt <- seq(0, 2, length.out = 200)
  y0 <- 1
  y1 <- deSolve::lsoda(y0, tt, function(...) list(target1(...)), numeric())
  y2 <- deSolve::lsoda(y0, -tt, function(...) list(target2(...)), numeric())
  expect_equal(y1[, 2], y2[, 2], tolerance = 1e-16)

  res1a <- dopri(y0, tt,  target1, numeric(0), return_statistics = TRUE)
  res2a <- dopri(y0, -tt, target2, numeric(0), return_statistics = TRUE)

  ## Not terribly accurate because of the nasty discontinuity
  expect_equal(res1a[, 2], y1[, 2], tolerance = 1e-4)
  ## As above, the dopri solutions should agree well
  expect_equal(res1a[, 2], res2a[, 2], tolerance = 1e-16)

  ## Add the critical time in:
  res1b <- dopri(y0, tt,  target1, numeric(0), tcrit = 1,
                 return_statistics = TRUE)
  res2b <- dopri(y0, -tt, target2, numeric(0), tcrit = -1,
                 return_statistics = TRUE)
  expect_equal(res1b[, 2], res2b[, 2], tolerance = 1e-16)

  s1a <- attr(res1a, "statistics")
  s2a <- attr(res2a, "statistics")
  s1b <- attr(res1b, "statistics")
  s2b <- attr(res2b, "statistics")
  expect_equal(s1a, s2a)
  expect_equal(s1b, s2b)
  expect_lt(s2b[["n_step"]], s2a[["n_step"]])
})

test_that("time validation", {
  expect_error(run_lorenz_dde(numeric(0)), "At least two times must be given")
  expect_error(run_lorenz_dde(0), "At least two times must be given")

  expect_error(run_lorenz_dde(c(0, 0)),
               "Beginning and end times are the same")
  expect_error(run_lorenz_dde(c(0, 2, 1)),
               "Times have inconsistent sign")
})

test_that("non-real input", {
  tt <- seq(0, 10, length.out = 101)
  cmp <- run_lorenz_dde(tt)
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  res <- dopri(y, tt, "lorenz", p, atol = 1e-7, rtol = 1e-7,
               parms_are_real = FALSE,
               dllname = "dde_lorenz2")
  expect_identical(res, cmp)
})

test_that("externalptr input", {
  tt <- seq(0, 10, length.out = 101)
  cmp <- run_lorenz_dde(tt)

  ptr <- .Call("lorenz_init",  c(10, 28, 8 / 3), PACKAGE = "dde_lorenz3")
  expect_is(ptr, "externalptr")
  y <- c(10, 1, 1)
  res <- dopri(y, tt, "lorenz", ptr, atol = 1e-7, rtol = 1e-7,
               parms_are_real = FALSE,
               dllname = "dde_lorenz3")
  expect_identical(res, cmp)
})

test_that("grow history", {
  tt <- seq(0, 10, length.out = 101)
  res <- run_lorenz_dde(tt, n_history = 5, grow_history = TRUE,
                        return_history = TRUE, return_statistics = TRUE)
  s <- attr(res, "statistics")
  h <- attr(res, "history")

  expect_equal(s[["n_accept"]], ncol(h))

  cmp <- run_lorenz_dde(tt, n_history = s[["n_step"]] * 2,
                        return_history = TRUE, return_statistics = TRUE)
  expect_equal(h, attr(cmp, "history"))
  expect_equal(s, attr(cmp, "statistics"))
})


test_that("initial derivative validation", {
  deriv <- function(t, y, p) {
    dy <- rep(1, length(y))
    if (length(p) > 0L) {
      dy[p] <- NA_real_
    }
    dy
  }

  tt <- c(0, 1)
  y0 <- rep(1, 10)

  expect_error(dopri(y0, tt, deriv, 1),
               "non-finite derivative at initial time for element 1")
  expect_error(dopri(y0, tt, deriv, 4:9),
               "non-finite derivative at initial time for element 4")
})


test_that("verbose mode prints trace", {
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
  names(y0) <- letters[1:3]

  expect_silent(ans1 <- dopri(y0, tt, lorenz, p))
  expect_output(ans2 <- dopri(y0, tt, lorenz, p, verbose = TRUE),
                ".step. t: .+, h: .+\n")
  expect_output(ans3 <- dopri(y0, tt, lorenz, p, verbose = VERBOSE_EVAL),
                ".eval. t: .+, h: .+\n")

  expect_identical(ans1, ans2)
  expect_identical(ans1, ans3)
})


test_that("verbose with callback", {
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
  names(y0) <- letters[1:3]

  steps <- list()
  evals <- list()
  callback <- function(t, h, y) {
    if (is.na(h)) {
      evals <<- c(evals, list(list(t = t, h = h, y = y)))
    } else {
      steps <<- c(steps, list(list(t = t, h = h, y = y)))
    }
  }

  ## Callback successfully returns all steps - we could have done
  ## something more interesting here, and typically this would involve
  ## printing not collection.
  ans1 <- dopri(y0, tt, lorenz, p, verbose = TRUE, callback = callback)
  n <- length(steps)
  expect_gt(n, 0)
  expect_equal(steps[[1]]$t, 0)
  expect_equal(steps[[1]]$y, unname(y0))
  expect_equal(steps[[n]]$t + steps[[n]]$h, 1)
  expect_equal(length(evals), 0)

  ans2 <- dopri(y0, tt, lorenz, p, verbose = VERBOSE_EVAL, callback = callback)
  expect_identical(steps[1:n], steps[-(1:n)])
  m <- length(evals)
  expect_gt(m, 0)
  expect_equal(evals[[1]]$t, 0)
  expect_equal(evals[[1]]$y, unname(y0))
  expect_true(all(is.na(vapply(evals, function(x) x$h, numeric(1)))))

  ## Time must have changed
  expect_gt(evals[[2]]$t, evals[[1]]$t)
  ## Output must have changed
  expect_true(all(evals[[2]]$y != evals[[1]]$y))
})


test_that("check verbose argument", {
  expect_identical(dopri_verbose(FALSE), VERBOSE_QUIET)
  expect_identical(dopri_verbose(0.0), VERBOSE_QUIET)
  expect_identical(dopri_verbose(0L), VERBOSE_QUIET)

  expect_identical(dopri_verbose(TRUE), VERBOSE_STEP)
  expect_identical(dopri_verbose(1.0), VERBOSE_STEP)
  expect_identical(dopri_verbose(1L), VERBOSE_STEP)

  expect_identical(dopri_verbose(2.0), VERBOSE_EVAL)
  expect_identical(dopri_verbose(2L), VERBOSE_EVAL)

  expect_error(dopri_verbose(-1), "Invalid value for verbose")
  expect_error(dopri_verbose(1.5), "verbose must be integer")
  expect_error(dopri_verbose(c(1, 2)), "verbose must be a scalar")
  expect_error(dopri_verbose(integer(0)), "verbose must be a scalar")
  expect_error(dopri_verbose(NULL), "verbose must be a scalar")
})


test_that("check callback argument", {
  expect_null(dopri_callback(NULL))

  f <- function(a, b, c) NULL
  res <- dopri_callback(f)
  expect_identical(f, res[[1]])
  expect_is(res[[2]], "environment")
  expect_equal(parent.env(res[[2]]), environment(f))

  expect_error(dopri_callback(TRUE),
               "Expected a function for 'callback'")
  expect_error(dopri_callback(function() 1),
               "Expected a function with 3 arguments for 'callback'")
})
