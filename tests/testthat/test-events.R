context("events")

test_that("events", {
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

  ## OK, this *looks* pretty good but we need to show that it is
  ## correct.  Once this is OK we have the first commit on this branch
  ## I think!

  ## This does seem to throw the tolerances out of whack quite badly,
  ## but then they're being computed on totally different quantities.
  m <- 2^findInterval(tt, events2$time + 1e-8,
                      rightmost.closed = TRUE)
  expect_equal(res2[, -1], cmp2[, -1] * m, tolerance = 1e-5)
})
