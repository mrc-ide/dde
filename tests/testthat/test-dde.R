context("dde")

test_that("dde", {
  skip_if_not_installed("deSolve")
  ## The pre-lag part of the integration
  tt1 <- seq(0, 14, length.out = 101)
  ## Post lag, but the first bit
  tt2 <- seq(0, 30, length.out = 301)
  ## The entire interesting phase
  tt3 <- seq(0, 200, length.out = 301)

  yy1 <- run_seir_deSolve(tt1)
  yy2 <- run_seir_deSolve(tt2)
  yy3 <- run_seir_deSolve(tt3)

  for (method in dopri_methods()) {
    ## Before the lag; this one is easy:
    expect_equal(run_seir_dde(tt1, method = method), yy1)

    expect_equal(run_seir_dde(tt2, method = method), yy2)

    ## Entire interesting region
    expect_equal(run_seir_dde(tt3, method = method), yy3, tolerance = 1e-6)
  }

  ## Confirm that different integrators were actually run here:
  y5 <- run_seir_dde(tt3, return_statistics = TRUE)
  y8 <- run_seir_dde(tt3, method = "dopri853", return_statistics = TRUE)
  expect_false(identical(y5, y8))

  ## The 853 stepper should take fewer steps (though in this case it's
  ## more ~10% more function evaluations because of more rejected
  ## steps).
  s5 <- attr(y5, "statistics")
  s8 <- attr(y8, "statistics")
  expect_gt(s5[["n_step"]], s8[["n_step"]])
  expect_lt(s5[["n_eval"]], s8[["n_eval"]])

  ## Run again with a critical time at the point the delay starts:
  y5_2 <- run_seir_dde(tt3, return_statistics = TRUE, tcrit = 14)
  y8_2 <- run_seir_dde(tt3, method = "dopri853", return_statistics = TRUE,
                       tcrit = 14)
  s5_2 <- attr(y5_2, "statistics")
  s8_2 <- attr(y8_2, "statistics")
  expect_gt(s5_2[["n_step"]], s8_2[["n_step"]])
  expect_gt(s5_2[["n_eval"]], s8_2[["n_eval"]])
  expect_lt(s8_2[["n_reject"]], s8[["n_reject"]])
})

test_that("output", {
  tt <- seq(0, 200, length.out = 301)
  p <- numeric(0)
  y0 <- c(1e7 - 1, 0, 1, 0)
  res1 <- dopri(y0, tt, "seir", p, n_history = 1000L,
                dllname = "dde_seir", return_history = FALSE)
  expect_equal(names(attributes(res1)), "dim")

  res2 <- dopri(y0, tt, "seir", p, n_history = 1000L,
                n_out = 1L, output = "seir_output",
                dllname = "dde_seir", return_history = FALSE,
                return_output_with_y = FALSE)

  output <- attr(res2, "output")
  expect_equal(dim(output), c(length(tt), 1L))
  expect_equal(drop(output), rowSums(res2[, -1]), tolerance = 1e-14)

  attr(res2, "output") <- NULL
  expect_identical(res1, res2)

  ## Corner case with the first output entry
  res3 <- dopri(y0, tt, "seir", p, n_history = 1000L,
                n_out = 1L, output = "seir_output",
                dllname = "dde_seir", return_history = FALSE,
                return_initial = FALSE, return_output_with_y = FALSE)
  m3 <- res3[]
  attr(m3, "output") <- NULL
  expect_identical(m3, res1[-1, ])
  out3 <- attr(res3, "output")
  expect_identical(out3, output[-1, , drop = FALSE])
})

test_that("R interface", {
  tt <- seq(0, 200, length.out = 301)
  p <- numeric(0)
  y0 <- c(1e7 - 1, 0, 1, 0)
  res1 <- dopri(y0, tt, "seir", p, n_history = 1000L,
                dllname = "dde_seir", return_history = FALSE)

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

    if (p == "nonnumeric") {
      tau <- p
    } else {
      tau <- t - lat_hum
      if (p == "as.integer" && abs(tau - as.integer(tau)) < 1e-6) {
        tau <- as.integer(tau)
      } else if (p == "nonscalar") {
        tau <- rep(tau, 2)
      }
    }
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

  for (method in dopri_methods()) {
    res2 <- dopri(y0, tt, seir, "one", n_history = 1000L,
                  return_history = FALSE, method = method)
    res3 <- dopri(y0, tt, seir, "idx", n_history = 1000L,
                  return_history = FALSE, method = method)
    res4 <- dopri(y0, tt, seir, "all", n_history = 1000L,
                  return_history = FALSE, method = method)

    expect_equal(res2, res1,
                 tolerance = if (method == "dopri5") 1e-8 else 1e-5)
    expect_identical(res3, res2)
    expect_identical(res4, res2)
  }

  ## Check some invalid input here too
  expect_error(dopri(y0, tt, seir, "nonnumeric", n_history = 1000L,
                     return_history = FALSE, method = method),
               "Expected a double")
  expect_error(dopri(y0, tt, seir, "nonscalar", n_history = 1000L,
                     return_history = FALSE, method = method),
               "Expected a scalar")
  res <- dopri(y0, tt, seir, "as.integer", n_history = 1000L,
               return_history = FALSE, method = method)
  expect_equal(res, res2)
})

test_that("R interface with output", {
  growth <- function(t, y, p) {
    tau <- t - 2.0
    c(y[[1L]], ylag(tau, 2L)) * p
  }
  output <- function(t, y, p) {
    tau <- t - 2.0
    ylag(tau, 1L)
  }
  tt <- seq(0, 20, length.out = 501)
  p <- 0.1
  y0 <- c(.1, .1)
  for (method in dopri_methods()) {
    if (method == "dopri5") {
      res <- dopri(y0, tt, growth, p,
                   n_out = 1, output = output,
                   n_history = 1000L, return_output_with_y = FALSE,
                   atol = 1e-8, rtol = 1e-8,
                   tcrit = 2, method = method)
    } else {
      expect_error(dopri(y0, tt, growth, p,
                         n_out = 1, output = output,
                         n_history = 1000L, return_output_with_y = FALSE,
                         tcrit = 2, method = method),
                   "step size vanished")
      ## The fix here is to include a number of multiples of the delay
      ## time that causes the discontinutites in the solution and all
      ## the bits that the integration relies on; this should be dealt
      ## with using some sort of discontinuty pruning approach I think.
      tcrit <- seq(2, 20, by = 2)
      res <- dopri(y0, tt, growth, p,
                   n_out = 1, output = output,
                   n_history = 1000L, return_output_with_y = FALSE,
                   atol = 1e-8, rtol = 1e-8,
                   tcrit = tcrit, method = method)
    }

    tol <- 1e-7
    i <- tt <= 2.0
    ## The first entry is easy:
    expect_equal(res[, 2], y0[1] * exp(p * tt), tolerance = tol)
    ## The second entry is in two parts, and only the first part is easy
    ## (for me) to compute:
    expect_equal(res[i, 3], y0[2] + y0[2] * p * tt[i], tolerance = tol)
    ## The output:
    out <- drop(attr(res, "output"))
    expect_equal(out[i], rep(y0[2], sum(i)))
    expect_equal(out[!i], y0[1] * exp(p * (tt[!i] - 2.0)), tolerance = tol)
  }
})

## There's a lot of really nasty bits here, so I'll build this test up
## in pieces.  Eventually we want to show that running a DDE with
## negative time gives the same result as an equivalent DDE with
## positive time.
test_that("delay, negative time", {
  growth0 <- function(t, y, p) {
    y * p
  }
  output0 <- function(t, y, p) {
    ylag(t - sign(p[1]) * 2)
  }

  set.seed(1)
  tt <- seq(0, 5, length.out = 501)
  r <- runif(2)
  y0 <- runif(2)
  ih <- seq_len(length(y0) * 5)

  real_fwd <- t(y0 * exp(outer(-r, -tt)))
  real_back <- t(y0 * exp(outer(r, tt)))
  expect_equal(real_fwd, real_back)

  ## First, the simplest case; running with no output
  res_fwd <- dopri(y0, tt, growth0, r, n_history = 1000L,
                   atol = 1e-8, rtol = 1e-8)
  res_back <- dopri(y0, -tt, growth0, -r, n_history = 1000L,
                   atol = 1e-8, rtol = 1e-8)

  expect_true(all.equal(res_fwd[, -1], real_fwd, check.attributes = FALSE))
  expect_true(all.equal(res_back[, -1], real_back, check.attributes = FALSE))
  expect_identical(attr(res_fwd, "history")[ih, ],
                   attr(res_back, "history")[ih, ])

  res_fwd <- dopri(y0, tt, growth0, r, output = output0, n_out = 2L,
                   n_history = 1000L, return_output_with_y = FALSE,
                   atol = 1e-8, rtol = 1e-8)
  res_back <- dopri(y0, -tt, growth0, -r, output = output0, n_out = 2L,
                    n_history = 1000L, return_output_with_y = FALSE,
                    atol = 1e-8, rtol = 1e-8)

  ## After this, the output agrees:
  real_output <- matrix(rep(y0, length(tt)), length(y0))
  i <- tt > 2
  real_output[, i] <- y0 * exp(outer(r, tt[i] - 2))
  real_output <- t(real_output)

  out_fwd <- attr(res_fwd, "output")
  out_back <- attr(res_back, "output")

  expect_true(all.equal(out_fwd, real_output))
  expect_equal(out_fwd, out_back)
})

test_that("failure to fetch history", {
  tt <- seq(0, 30, length.out = 301)
  expect_error(run_seir_dde(tt, n_history = 2L),
               "Integration failure: did not find time")
})

test_that("ylag_vec_int", {
  times <- seq(0, 20, length.out = 101)
  p <- numeric(0)
  y0 <- c(1e7 - 1, 0, 1, 0)
  for (method in dopri_methods()) {
    cmp <- dopri(y0, times, "seir", p, n_history = 1000L, method = method,
                 dllname = "dde_seir")
    res <- dopri(y0, times, "seir", p, n_history = 1000L, method = method,
                 dllname = "dde_seir_int")
    expect_equal(cmp, res, tolerance = 1e-8)
  }
})

test_that("Zero lag time", {
  ## Oops; this is not good; delay is too short
  growth <- function(t, y, p) {
    c(y[[1L]], ylag(t, 2L)) * p
  }
  tt <- seq(0, 20, length.out = 501)
  p <- 0.1
  y0 <- c(.1, .1)
  expect_error(dopri(y0, tt, growth, p, n_history = 1000L),
               "did not find time in history")
})

## OK, most of the fundamentals are here now, but need to get the
## argument handling correct (I'm passing a NULL through where there
## should have been a TRUE/FALSE).  Once this is basically working,
## I'll get this going for the difeq version too.
test_that("restart", {
  tt <- seq(0, 200, length.out = 101)
  tt1 <- tt[tt < 80]
  tt2 <- tt[tt >= tt1[length(tt1)]]

  cmp <- run_seir_dde(tt, return_minimal = TRUE)
  cmp1 <- run_seir_dde(tt1, return_minimal = TRUE)
  cmp2 <- cmp[, tt[-1] >= tt1[length(tt1)]]

  res1 <- run_seir_dde(tt1, restartable = TRUE, return_minimal = TRUE)
  expect_is(attr(res1, "ptr"), "externalptr")

  expect_equal(res1, cmp1, check.attributes = FALSE)

  ## Then start work on the continuation
  res2 <- dopri_continue(res1, tt2, return_initial = TRUE)

  ## NOTE: testthat will report incorrectly here on failure because of
  ## attributes (at least of 1.0.2)
  expect_equal(res2, cmp2, check.attributes = FALSE)

  expect_error(dopri_continue(res1, tt2, return_initial = TRUE),
               "Incorrect initial time on integration restart")
})

test_that("restart and copy", {
  tt <- seq(0, 200, length.out = 101)
  tt1 <- tt[tt < 80]
  tt2 <- tt[tt >= tt1[length(tt1)]]

  cmp <- run_seir_dde(tt, return_minimal = TRUE)
  cmp1 <- run_seir_dde(tt1, return_minimal = TRUE)
  cmp2 <- cmp[, tt[-1] >= tt1[length(tt1)]]

  res1 <- run_seir_dde(tt1, restartable = TRUE, return_minimal = TRUE)
  expect_is(attr(res1, "ptr"), "externalptr")

  res2 <- dopri_continue(res1, tt2, return_initial = TRUE, copy = TRUE)
  res3 <- dopri_continue(res1, tt2, return_initial = TRUE, copy = TRUE)
  expect_equal(res2, cmp2, check.attributes = FALSE)
  expect_equal(res2, res3)

  ## Continue off and things should work OK here
  res4 <- dopri_continue(res1, tt2, return_initial = TRUE, copy = FALSE)
  expect_equal(res4, cmp2, check.attributes = FALSE)

  ## But we've modified the pointer so this will no longer work:
  expect_error(dopri_continue(res1, tt2, return_initial = TRUE),
               "Incorrect initial time on integration restart")
})

test_that("restart errors", {
  tt <- seq(0, 200, length.out = 101)
  tt1 <- tt[tt < 80]
  tt2 <- tt[tt >= tt1[length(tt1)]]

  cmp <- run_seir_dde(tt, return_minimal = TRUE)
  cmp1 <- run_seir_dde(tt1, return_minimal = TRUE)
  cmp2 <- cmp[, tt[-1] >= tt1[length(tt1)]]

  res1 <- run_seir_dde(tt1, restartable = TRUE, return_minimal = TRUE)

  expect_error(dopri_continue(res1, tt2[[1]], copy = TRUE),
               "At least two times must be given")
  expect_error(dopri_continue(res1, rev(tt1), copy = TRUE),
               "Incorrect sign for the times")

  y1 <- res1[, ncol(res1)]
  expect_error(dopri_continue(res1, tt2, y1[-1], copy = TRUE),
               "Incorrect size 'y' on integration restart")
  expect_error(dopri_continue(res1, tt2, rep(y1, 2), copy = TRUE),
               "Incorrect size 'y' on integration restart")

  ptr <- attr(res1, "ptr")
  ## Make a NULL pointer:
  attr(res1, "ptr") <- make_null_pointer(ptr)
  expect_error(dopri_continue(res1, tt2, copy = TRUE),
               "pointer has been freed")
  ## Delete the pointer:
  attr(res1, "ptr") <- NULL
  expect_error(dopri_continue(res1, tt2, copy = TRUE),
               "Expected an external pointer")
  ## Something stupid as a pointer:
  attr(res1, "ptr") <- res1
  expect_error(dopri_continue(res1, tt2, copy = TRUE),
               "Expected an external pointer")
})

test_that("change y on restart", {
  growth <- function(t, y, p) {
    y * p
  }
  output <- function(t, y, p) {
    ylag(t - 2.0)
  }

  tt <- seq(0, 10, length.out = 101)
  tc <- 4
  i1 <- tt <= tc
  i2 <- tt >= tc
  tt1 <- tt[i1]
  tt2 <- tt[i2]

  y0 <- 1
  r <- -0.5
  j <- seq_along(y0)

  ## TODO: error is not very clear when output = NULL but n_out > 0
  ##
  ## TODO: error is not very clear when n_history = 0
  ##
  ## TODO: error if tcrit[1] == t[1]

  res <- dopri(y0, tt, growth, r, n_out = length(y0), output = output,
               n_history = 1000, return_history = FALSE,
               return_time = FALSE, return_by_column = FALSE,
               tcrit = tt2[[1]], atol = 1e-8, rtol = 1e-8)

  res1 <- dopri(y0, tt1, growth, r, n_out = length(y0), output = output,
                n_history = 1000, return_history = FALSE,
                return_time = FALSE, return_by_column = FALSE,
                restartable = TRUE)

  y1 <- res1[j, ncol(res1)]
  res2 <- dopri_continue(res1, tt2, y1, copy = TRUE)

  ## Check that options have been preserved:
  expect_equal(nrow(res2), nrow(res1))
  expect_equal(ncol(res2), length(tt2))
  expect_null(attr(res2, "output"))

  ## TODO: not sure why this is very slightly off, but it's possibly
  ## due to a change in step size?  I don't think that this should be
  ## the case though and would like to see this shrink.
  expect_equal(res2[j, ], res[j, i2], tolerance = 5e-7)
  expect_equal(res2[-j, ], res[-j, i2], tolerance = 5e-7)

  ## Change y on re-entry:
  y2 <- y1 * 2
  res3 <- dopri_continue(res1, tt2, y2, copy = TRUE, tcrit = tt2[[1]] + 2)

  ## TODO: this is not terrific accuracy either:
  expect_equal(res3[j, ] / res2[j, ], rep(2, length(tt2)), tolerance = 1e-5)
  k <- tt2 - tt2[[1]] < 2.0
  expect_equal(res3[-j, k] / res2[-j, k], rep(1, sum(k)))
  expect_equal(res3[-j, !k] / res2[-j, !k], rep(2, sum(!k)), tolerance = 1e-5)
})

## This indicates to me that we might still have a problem with the
## solver, but it might also be that the high order difficulty from
## the echo of the delay is causing trouble.  If you push out tcrit
## values as multiples of 14 and drop tolerance to 1e-5 you can get as
## far as t = 56 before it gives up.
test_that("stiffness detection", {
  tt <- seq(0, 200, length.out = 301)
  yy <- run_seir_dde(tt)
  expect_equal(run_seir_dde(tt, stiff_check = 1), yy)
  expect_error(run_seir_dde(tt, method = "dopri853", stiff_check = 1),
               "problem became stiff")
})
