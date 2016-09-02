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

test_that("integer time input", {
  tt <- seq(0, 10, by=1)
  expect_identical(run_lorenz_dde(as.integer(tt)),
                   run_lorenz_dde(tt))
})

test_that("zero time difference", {
  tt <- c(0, 0, 1, 1, 2, 2, 2, 3, 3, 3)
  res1 <- run_lorenz_dde(tt, return_initial=TRUE)
  res2 <- run_lorenz_dde(unique(tt), return_initial=TRUE)
  expect_identical(res1, res2[, tapply(tt, tt)])
})

test_that("ode, 873 stepper", {
  tt <- seq(0, 1, length.out = 200)
  m5 <- run_lorenz_dde(tt)
  m8 <- run_lorenz_dde(tt, method = "dopri853")
  expect_equal(m5, m8, tolerance=1e-7)
  expect_false(identical(m5, m8))
})

test_that("dense output", {
  tt <- seq(0, 1, length.out = 200)
  m1 <- run_lorenz_deSolve(tt)
  dimnames(m1) <- NULL

  for (method in c(dopri_methods())) {
    m2 <- run_lorenz_dde(c(0, max(tt)), n_history = 1000L, method = method)
    ## single row of output:
    expect_equal(dim(m2), c(3, 1))
    h2 <- attr(m2, "history")
    expect_equal(nrow(h2),
                 if (method == "dopri5") 17L else 26L) # c(5, 8) * 3 + 2
    expect_identical(attr(h2, "n"), 3L)

    m3 <- dopri_interpolate(h2, tt)
    m4 <- run_lorenz_dde(tt, method = method)

    expect_identical(m3[-1,], t(m4))
    expect_equal(m3, m1, tolerance=1e-6)

    ## Check column output:
    m5 <- run_lorenz_dde(tt, n_history = 1000L, by_column = TRUE, method = method)
    expect_identical(attr(m5, "history"), attr(m2, "history"))
    attr(m5, "history") <- NULL
    expect_identical(m5, t(m4))

    expect_error(dopri_interpolate(h2, tt - 1),
                 "Time falls outside of range of known history")
    expect_error(dopri_interpolate(h2, tt + 1),
                 "Time falls outside of range of known history")

    expect_error(dopri_interpolate(h2[-1,], tt),
                 "Corrupt history object")
    h2 <- h2[-1,]
    attr(h2, "n") <- 3L
    expect_error(dopri_interpolate(h2, tt),
                 "Corrupt history object: incorrect number of rows")
  }
})

test_that("output", {
  tt <- seq(0, 1, length.out = 200)

  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)

  res1 <- dopri(y0, tt, "lorenz", p, dllname = "lorenz")
  res2 <- dopri(y0, tt, "lorenz", p, dllname = "lorenz",
                n_out = 2L, output = "lorenz_output")

  expect_equal(names(attributes(res1)), "dim")
  output <- attr(res2, "output")
  expect_equal(dim(output), c(2L, ncol(res1)))
  expect_equal(output[1L, ], pmin(res1[1,], res1[2,], res1[3,]))
  expect_equal(output[2L, ], pmax(res1[1,], res1[2,], res1[3,]))

  attr(res2, "output") <- NULL
  expect_identical(res1, res2)

  res3 <- dopri(y0, tt, "lorenz", p, dllname = "lorenz",
                n_out = 2L, output = "lorenz_output", by_column = TRUE)
  expect_equal(attr(res3, "output"), t(output))
  attr(res3, "output") <- NULL
  expect_identical(res3, t(res2))
})

test_that("keep initial", {
  tt <- seq(0, 1, length.out = 200)

  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)

  res1 <- dopri(y0, tt, "lorenz", p, dllname = "lorenz",
                return_initial = TRUE)
  res2 <- dopri(y0, tt, "lorenz", p, dllname = "lorenz",
                return_initial = FALSE)
  expect_equal(ncol(res1), length(tt))
  expect_identical(res1[, 1], y0)
  expect_identical(res1[, -1], res2)

  res3 <- dopri(y0, tt, "lorenz", p, dllname = "lorenz",
                n_out = 2L, output = "lorenz_output", return_initial = TRUE)
  res4 <- dopri(y0, tt, "lorenz", p, dllname = "lorenz",
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
  res1 <- dopri(y0, tt, "lorenz", p, dllname = "lorenz")
  res2 <- dopri(y0, tt, lorenz, p)
  expect_identical(res1, res2)

  res3 <- dopri(y0, tt, "lorenz", p, dllname = "lorenz",
                n_out = 2L, output = "lorenz_output")
  res4 <- dopri(y0, tt, lorenz, p, n_out = 2L, output = lorenz_output)
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
  res1 <- dopri(1, tt, target, numeric(0), return_statistics = TRUE)
  res2 <- dopri(1, tt, target, numeric(0), tcrit=1, return_statistics = TRUE)

  s1 <- attr(res1, "statistics")
  s2 <- attr(res2, "statistics")
  expect_lt(s2[["n_step"]], s1[["n_step"]])

  expect_is(attr(res1, "step_size"), "numeric")
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

  expect_null(dimnames(dopri(y0, tt, lorenz, p)))
  expect_equal(dimnames(dopri(y0, tt, lorenz, p, ynames=nms)), cmp)
  expect_equal(dimnames(dopri(setNames(y0, nms), tt, lorenz, p)), cmp)
  expect_null(dimnames(dopri(setNames(y0, nms), tt, lorenz, p, ynames=FALSE)))

  expect_error(dopri(y0, tt, lorenz, p, ynames=nms[1]),
               "ynames must be the same length as y")
  expect_error(dopri(y0, tt, lorenz, p, ynames=1),
               "Invalid value for ynames")

  ## Similar for output names:
  onms <- LETTERS[1:2]
  ocmp <- list(onms, NULL)
  expect_null(dimnames(attr(dopri(y0, tt, lorenz, p, n_out = 2L,
                                  output = lorenz_output), "output")))
  expect_null(dimnames(attr(dopri(y0, tt, lorenz, p, n_out = 2L,
                                  output = lorenz_output, outnames = NULL),
                            "output")))
  expect_equal(dimnames(attr(dopri(y0, tt, lorenz, p, n_out = 2L,
                                   output = lorenz_output, outnames = onms),
                             "output")), ocmp)
  expect_error(dopri(y0, tt, lorenz, p, n_out = 2L,
                     output = lorenz_output, outnames = nms),
               "outnames must have length n_out")
  expect_error(dopri(y0, tt, lorenz, p, n_out = 2L,
                     output = lorenz_output, outnames = 1),
               "Invalid value for outnames")

  ## Check both together:
  res <- dopri(y0, tt, lorenz, p, n_out = 2L,
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

  expect_equal(nrow(dopri(y0, tt, lorenz, p)), 3)
  res <- dopri(y0, tt, lorenz, p, return_time = TRUE)
  expect_equal(nrow(res), 4)
  expect_equal(res[1, ], tt[-1L])

  ## include the first time:
  expect_equal(
    dopri(y0, tt, lorenz, p, return_time = TRUE, return_initial = TRUE)[1, ],
    tt)

  ## with names:
  nms <- letters[1:3]
  res <- dopri(y0, tt, lorenz, p, return_time = TRUE, ynames = nms)
  expect_equal(rownames(res), c("time", nms))

  res <- dopri(y0, tt, lorenz, p, return_time = TRUE, ynames = nms,
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

  res <- dopri(y0, tt, lorenz, p, deSolve_compatible = TRUE)
  expect_equal(dim(res), c(length(tt), 4))
  expect_equal(dimnames(res), list(NULL, c("time", 1:3)))
  expect_equal(res[, 1], tt)
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

test_that("Native Symbol interface", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 10, length.out=101)
  func <- getNativeSymbolInfo("lorenz", PACKAGE="lorenz")
  y <- dopri(y, times, func, p)
  expect_equal(y, run_lorenz_dde(times, tol = 1e-6))
})

test_that("NULL pointer safety", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 10, length.out=101)

  func <- getNativeSymbolInfo("lorenz", PACKAGE="lorenz")
  func <- unserialize(serialize(func, NULL))

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
    dopri(y0, tt, growth, r, n_history = 100, deSolve_compatible=TRUE)))
  expect_true(has_history(
    dopri(y0, tt, growth, r, n_history = 100, by_column=TRUE)))
  expect_true(has_history(
    dopri(y0, tt, growth, r, n_history = 100, return_initial=TRUE)))
  expect_true(has_history(
    dopri(y0, tt, growth, r, n_history = 100, return_time=TRUE)))

  expect_false(has_statistics(dopri(y0, tt, growth, r)))
  expect_true(has_statistics(dopri(y0, tt, growth, r,
                                   return_statistics = TRUE)))

  expect_true(has_statistics(
    dopri(y0, tt, growth, r, return_statistics = TRUE,
          deSolve_compatible=TRUE)))
  expect_true(has_statistics(
    dopri(y0, tt, growth, r, return_statistics = TRUE,
          by_column=TRUE)))
  expect_true(has_statistics(
    dopri(y0, tt, growth, r, return_statistics = TRUE,
          return_initial=TRUE)))
  expect_true(has_statistics(
    dopri(y0, tt, growth, r, return_statistics = TRUE,
          return_time=TRUE)))
})

test_that("negative time", {
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
  expect_equal(cmp[, -1], cmp_fwd[, -1], tolerance=1e-16)

  ## and dde:
  res <- dopri(y0, tt, growth, r, deSolve_compatible = TRUE,
               n_history = 100)
  res_fwd <- dopri(y0, -tt, growth, -r, deSolve_compatible = TRUE,
                   n_history = 100)

  expect_equal(res[, -1], res_fwd[, -1], tolerance=1e-16)

  ## Interestingly, we do stop at the same points here, so the error
  ## calculation and underlying stepping is probably OK.
  t_idx <- length(y0) * 5 + 1
  h1 <- attr(res, "history")
  h2 <- attr(res_fwd, "history")
  expect_equal(-h1[t_idx, ], h2[t_idx, ], tolerance=1e-16)
  expect_equal(-h1[t_idx + 1, ], h2[t_idx + 1, ], tolerance=1e-16)
  expect_equal(h1[1:20,], h2[1:20,], tolerance=1e-16)
})

test_that("negative time with tcrit", {
  ## TODO: To be a really good test, this would need to test the case
  ## with more than one critical time, as there are a couple of places
  ## where the critical time gets advanced.
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
  expect_equal(y1[, 2], y2[, 2], tolerance=1e-16)

  res1a <- dopri(y0, tt,  target1, numeric(0),
                 return_statistics = TRUE, return_initial=TRUE)
  res2a <- dopri(y0, -tt, target2, numeric(0),
                 return_statistics = TRUE, return_initial=TRUE)

  ## Not terribly accurate because of the nasty discontinuity
  expect_equal(res1a[1, ], y1[, 2], tolerance=1e-4)
  ## As above, the dopri solutions should agree well
  expect_equal(res1a[1, ], res2a[1, ], tolerance=1e-16)

  ## Add the critical time in:
  res1b <- dopri(y0, tt,  target1, numeric(0), tcrit = 1,
                 return_statistics = TRUE, return_initial = TRUE)
  res2b <- dopri(y0, -tt, target2, numeric(0), tcrit = -1,
                 return_statistics = TRUE, return_initial = TRUE)
  expect_equal(res1b[1, ], res2b[1, ], tolerance=1e-16)

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
