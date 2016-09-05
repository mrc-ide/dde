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
  expect_equal(run_seir_dde(tt1),
               t(yy1[-1, ]))
  expect_equal(run_seir_dde(tt1, method = "dopri853"),
               t(yy1[-1, ]))

  ## Post lag
  expect_equal(run_seir_dde(tt2),
               t(yy2[-1, ]))
  expect_equal(run_seir_dde(tt2, method = "dopri853"),
               t(yy2[-1, ]))

  ## Entire interesting
  expect_equal(run_seir_dde(tt3),
               t(yy3[-1, ]), tolerance = 1e-6)
  expect_equal(run_seir_dde(tt3, method="dopri853"),
               t(yy3[-1, ]), tolerance = 1e-6)

  ## Confirm that different integrators were actually run here:
  y5 <- run_seir_dde(tt3, return_statistics=TRUE)
  y8 <- run_seir_dde(tt3, method="dopri853", return_statistics=TRUE)
  expect_false(identical(y5, y8))

  ## The 853 stepper should take fewer steps (though in this case it's
  ## more ~10% more function evaluations because of more rejected
  ## steps).
  s5 <- attr(y5, "statistics")
  s8 <- attr(y8, "statistics")
  expect_gt(s5[["n_step"]], s8[["n_step"]])
  expect_lt(s5[["n_eval"]], s8[["n_eval"]])

  ## Run again with a critical time at the point the delay starts:
  y5_2 <- run_seir_dde(tt3, return_statistics=TRUE, tcrit = 14)
  y8_2 <- run_seir_dde(tt3, method="dopri853", return_statistics=TRUE, tcrit = 14)
  s5_2 <- attr(y5_2, "statistics")
  s8_2 <- attr(y8_2, "statistics")
  expect_gt(s5_2[["n_step"]], s8_2[["n_step"]])
  expect_gt(s5_2[["n_eval"]], s8_2[["n_eval"]])
  expect_lt(s8_2[["n_reject"]], s8[["n_reject"]])
})

test_that("output", {
  tt <- seq(0, 200, length.out=301)
  p <- numeric(0)
  y0 <- c(1e7 - 1, 0, 1, 0)
  res1 <- dopri(y0, tt, "seir", p, n_history = 1000L,
                dllname = "seir", return_history = FALSE)
  expect_equal(names(attributes(res1)), "dim")

  res2 <- dopri(y0, tt, "seir", p, n_history = 1000L,
                n_out = 1L, output = "seir_output",
                dllname = "seir", return_history = FALSE)

  output <- attr(res2, "output")
  expect_equal(dim(output), c(1L, ncol(res1)))
  expect_equal(drop(output), colSums(res2), tolerance = 1e-14)

  attr(res2, "output") <- NULL
  expect_identical(res1, res2)

  ## Corner case with the first output entry
  res3 <- dopri(y0, tt, "seir", p, n_history = 1000L,
                n_out = 1L, output = "seir_output",
                dllname = "seir", return_history = FALSE,
                return_initial = TRUE)
  expect_identical(res3[, -1], res1[])
  out3 <- attr(res3, "output")
  expect_identical(out3[, -1, drop=FALSE], output)
  expect_identical(out3[, 1], output[, 1])
})

test_that("R interface", {
  tt <- seq(0, 200, length.out=301)
  p <- numeric(0)
  y0 <- c(1e7 - 1, 0, 1, 0)
  res1 <- dopri(y0, tt, "seir", p, n_history = 1000L,
                dllname = "seir", return_history = FALSE)

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

  res2 <- dopri(y0, tt, seir, "one", n_history = 1000L, return_history = FALSE)
  res3 <- dopri(y0, tt, seir, "idx", n_history = 1000L, return_history = FALSE)
  res4 <- dopri(y0, tt, seir, "all", n_history = 1000L, return_history = FALSE)

  expect_equal(res2, res1, tolerance = 1e-14)
  expect_identical(res3, res2)
  expect_identical(res4, res2)
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
  tt <- seq(0, 20, length.out=501)
  p <- 0.1
  y0 <- c(.1, .1)
  method <- dopri_methods()[[2]]
  for (method in dopri_methods()) {
    if (method == "dopri5") {
      res <- dopri(y0, tt, growth, p,
                   n_out = 1, output = output,
                   n_history = 1000L, return_initial = TRUE,
                   atol = 1e-8, rtol = 1e-8,
                   tcrit = 2, method = method)
    } else {
      expect_error(dopri(y0, tt, growth, p,
                         n_out = 1, output = output,
                         n_history = 1000L, return_initial = TRUE,
                         tcrit = 2, method = method),
                   "step size vanished")
      ## The fix here is to include a number of multiples of the delay
      ## time that causes the discontinutites in the solution and all
      ## the bits that the integration relies on; this should be dealt
      ## with using some sort of discontinuty pruning approach I think.
      tcrit <- seq(2, 20, by=2)
      res <- dopri(y0, tt, growth, p,
                   n_out = 1, output = output,
                   n_history = 1000L, return_initial = TRUE,
                   atol = 1e-8, rtol = 1e-8,
                   tcrit = tcrit, method = method)
    }

    tol <- 1e-7
    i <- tt <= 2.0
    ## The first entry is easy:
    expect_equal(res[1, ], y0[1] * exp(p * tt), tolerance=tol)
    ## The second entry is in two parts, and only the first part is easy
    ## (for me) to compute:
    expect_equal(res[2, i], y0[2] + y0[2] * p * tt[i], tolerance=tol)
    ## The output:
    out <- drop(attr(res, "output"))
    expect_equal(out[i], rep(y0[2], sum(i)))
    expect_equal(out[!i], y0[1] * exp(p * (tt[!i] - 2.0)), tolerance=tol)
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
  tt <- seq(0, 5, length.out=501)
  r <- runif(2)
  y0 <- runif(2)
  ih <- seq_len(length(y0) * 5)

  real_fwd <- t(y0 * exp(outer(-r, -tt)))
  real_back <- t(y0 * exp(outer(r, tt)))
  expect_equal(real_fwd, real_back)

  ## First, the simplest case; running with no output
  res_fwd <- dopri(y0, tt, growth0, r,
                   n_history = 1000L, return_initial=TRUE,
                   atol = 1e-8, rtol = 1e-8)
  res_back <- dopri(y0, -tt, growth0, -r,
                   n_history = 1000L, return_initial=TRUE,
                   atol = 1e-8, rtol = 1e-8)

  expect_true(all.equal(t(res_fwd[]), real_fwd, check.attributes=FALSE))
  expect_true(all.equal(t(res_back[]), real_back, check.attributes=FALSE))
  expect_identical(attr(res_fwd, "history")[ih, ],
                   attr(res_back, "history")[ih, ])

  res_fwd <- dopri(y0, tt, growth0, r,
                   output = output0, n_out = 2L,
                   n_history = 1000L, return_initial=TRUE,
                   atol = 1e-8, rtol = 1e-8)
  res_back <- dopri(y0, -tt, growth0, -r,
                    output = output0, n_out = 2L,
                    n_history = 1000L, return_initial=TRUE,
                    atol = 1e-8, rtol = 1e-8)

  ## After this, the output agrees:
  real_output <- matrix(rep(y0, length(tt)), length(y0))
  i <- tt > 2
  real_output[, i] <- y0 * exp(outer(r, tt[i] - 2))

  out_fwd <- attr(res_fwd, "output")
  out_back <- attr(res_back, "output")

  expect_true(all.equal(out_fwd, real_output))
  expect_equal(out_fwd, out_back)
})

## Next, try a restart; we'll run a system with some history and save
## everything, then modify the system and do a restart.
