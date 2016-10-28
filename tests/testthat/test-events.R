context("events")

test_that("change variables", {
  target_r <- function(t, y, p) {
    -y * p
  }
  event_r <- function(t, y, p) {
    y * 2
  }
  target_c <- "exponential"
  event_c <- "double_variables"
  dllname <- "growth"

  set.seed(1)
  y0 <- runif(5)
  p <- rep(1, length(y0))
  tol <- 1e-8
  tt <- seq(0, 5, length.out = 51) # 5001

  cmp1_r <- dopri(y0, tt, target_r, p, atol = tol, rtol = tol)
  cmp2_r <- dopri(y0, tt, target_r, p, tcrit = 1:5, atol = tol, rtol = tol)
  cmp1_c <- dopri(y0, tt, target_c, p, dllname = dllname,
                  atol = tol, rtol = tol)
  cmp2_c <- dopri(y0, tt, target_c, p, dllname = dllname, tcrit = 1:5,
                  atol = tol, rtol = tol)

  expect_equal(cmp1_r, cmp2_r, tolerance = 1e-6)
  expect_equal(cmp1_r, cmp1_c, tolerance = 1e-12)
  expect_equal(cmp2_r, cmp2_c, tolerance = 1e-12)

  events1_r <- list(time = 1:5, event = function(t, y, p) y)
  events2_r <- list(time = 1:5, event = event_r)
  events1_c <- list(time = 1:5, event = "identity")
  events2_c <- list(time = 1:5, event = event_c)

  ## First, the dummy event:
  res1_r <- dopri(y0, tt, target_r, p, events = events1_r,
                  atol = tol, rtol = tol)
  expect_equal(res1_r, cmp2_r, tolerance = 1e-12)

  ## Then the real event where we double the output variable:
  res2_r <- dopri(y0, tt, target_r, p, events = events2_r,
                  atol = tol, rtol = tol)

  ## This does seem to throw the tolerances out of whack quite badly,
  ## but then they're being computed on totally different quantities.
  m <- 2^findInterval(tt, events2_r$time + 1e-8,
                      rightmost.closed = TRUE)
  expect_equal(res2_r[, -1], cmp2_r[, -1] * m, tolerance = 1e-5)

  ## Confirm that the C cases agree fairly well:
  res1_c <- dopri(y0, tt, target_c, p, events = events1_c, dllname = dllname,
                  atol = tol, rtol = tol)
  expect_equal(res1_c, res1_r, tolerance = 1e-12)
  res2_c <- dopri(y0, tt, target_c, p, events = events2_c, dllname = dllname,
                  atol = tol, rtol = tol)
  expect_equal(res2_c, res2_r, tolerance = 1e-12)
})

test_that("change parameters", {
  target_r <- function(t, y, p) {
    p
  }
  event_r <- function(t, y, p) {
    attr(y, "parms") <- p * 2
    y
  }

  target_c <- "linear"
  event_c <- "double_parameters"
  dllname <- "growth"

  set.seed(1)
  y0 <- runif(5)
  p <- runif(5)
  p0 <- p
  tol <- 1e-8
  tt <- seq(0, 5, length.out = 51) # 5001

  cmp1_r <- dopri(y0, tt, target_r, p, atol = tol, rtol = tol)
  cmp2_r <- dopri(y0, tt, target_r, p, tcrit = 1:5, atol = tol, rtol = tol)
  cmp1_c <- dopri(y0, tt, target_c, p, dllname = dllname,
                  atol = tol, rtol = tol)
  cmp2_c <- dopri(y0, tt, target_c, p, dllname = dllname, tcrit = 1:5,
                  atol = tol, rtol = tol)

  expect_equal(cmp1_r, cmp2_r, tolerance = 1e-6)
  expect_equal(cmp1_r, cmp1_c, tolerance = 1e-12)
  expect_equal(cmp2_r, cmp2_c, tolerance = 1e-12)

  events1_r <- list(time = 1:5, event = function(t, y, p) y)
  events2_r <- list(time = 1:5, event = event_r)
  events1_c <- list(time = 1:5, event = "identity")
  events2_c <- list(time = 1:5, event = event_c)

  ## First, the dummy event (duplicates the test above really)
  res1_r <- dopri(y0, tt, target_r, p, events = events1_r,
                  atol = tol, rtol = tol)
  expect_equal(res1_r, cmp2_r, tolerance = 1e-12)

  ## Then the real event where we double the parameters
  res2_r <- dopri(y0, tt, target_r, p, events = events2_r,
                  atol = tol, rtol = tol)
  expect_identical(p, p0) # unchanged -- this is important!

  cmp <- matrix(NA, length(tt), length(y0))
  cmp[1, ] <- y0
  ## Then start the process of computing the trajectory.  This is
  ## actually quite a faff; surprisingly more than the above because I
  ## don't know how to do this sort of thing very effectively.
  times <- c(tt[1], events2_r$time)
  for (i in seq_along(events2_r$time)) {
    j <- which(tt >= times[i] & tt <= times[i + 1])
    cmp[j, ] <- t(cmp[j[1], ] + outer(p * 2 ^ (i - 1), tt[j] - tt[j][1]))
  }

  expect_equal(res2_r[, -1], cmp, tolerance = 1e-5)

  ## Confirm that the C cases agree fairly well:
  res1_c <- dopri(y0, tt, target_c, p, events = events1_c, dllname = dllname,
                  atol = tol, rtol = tol)
  expect_equal(res1_c, res1_r, tolerance = 1e-12)
  res2_c <- dopri(y0, tt, target_c, p, events = events2_c, dllname = dllname,
                  atol = tol, rtol = tol)
  expect_equal(res2_c, res2_r, tolerance = 1e-12)
})
