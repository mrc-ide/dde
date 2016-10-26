context("events")

test_that("change variables (R)", {
  target <- function(t, y, p) {
    -y
  }
  myevent <- function(t, y, p) {
    y * 2
  }

  set.seed(1)
  y0 <- runif(5)
  p <- NULL
  tol <- 1e-8
  tt <- seq(0, 5, length.out = 51) # 5001
  cmp1 <- dopri(y0, tt, target, p, atol = tol, rtol = tol)
  cmp2 <- dopri(y0, tt, target, p, tcrit = 1:5, atol = tol, rtol = tol)

  expect_equal(cmp1, cmp2, tolerance = 1e-6)

  events1 <- list(time = 1:5, event = function(t, y, p) y)
  res1 <- dopri(y0, tt, target, p, events = events1, atol = tol, rtol = tol)
  expect_equal(res1, cmp2, tolerance = 1e-12)

  events2 <- list(time = 1:5, event = myevent)
  res2 <- dopri(y0, tt, target, p, events = events2, atol = tol, rtol = tol)

  ## This does seem to throw the tolerances out of whack quite badly,
  ## but then they're being computed on totally different quantities.
  m <- 2^findInterval(tt, events2$time + 1e-8,
                      rightmost.closed = TRUE)
  expect_equal(res2[, -1], cmp2[, -1] * m, tolerance = 1e-5)
})

test_that("change parameters (R)", {
  target <- function(t, y, p) {
    p
  }
  myevent <- function(t, y, p) {
    attr(y, "parms") <- p * 2
    y
  }

  set.seed(1)
  y0 <- runif(5)
  p <- runif(5)
  p0 <- p
  tol <- 1e-8
  tt <- seq(0, 5, length.out = 51) # 5001

  ## This one is a bit harder.
  events1 <- list(time = 1:5, event = function(t, y, p) y)
  events2 <- list(time = 1:5, event = myevent)

  cmp1 <- dopri(y0, tt, target, p, atol = tol, rtol = tol)
  cmp2 <- dopri(y0, tt, target, p, tcrit = 1:5, atol = tol, rtol = tol)

  expect_equal(cmp1, cmp2, tolerance = 1e-6)

  res1 <- dopri(y0, tt, target, p, events = events1, atol = tol, rtol = tol)
  expect_equal(res1, cmp2, tolerance = 1e-12)

  res2 <- dopri(y0, tt, target, p, events = events2, atol = tol, rtol = tol)
  expect_identical(p, p0) # unchanged -- this is important!

  cmp <- matrix(NA, length(tt), length(y0))
  cmp[1, ] <- y0
  ## Then start the process of computing the trajectory.  This is
  ## actually quite a faff; surprisingly more than the above because I
  ## don't know how to do this sort of thing very effectively.
  times <- c(tt[1], events2$time)
  for (i in seq_along(events2$time)) {
    j <- which(tt >= times[i] & tt <= times[i + 1])
    cmp[j, ] <- t(cmp[j[1], ] + outer(p * 2 ^ (i - 1), tt[j] - tt[j][1]))
  }

  expect_equal(res2[, -1], cmp, tolerance = 1e-5)
})
