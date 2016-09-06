context("discrete")

## Solve some basic discrete time equations
test_that("increase", {
  rhs <- function(i, t, y, p) {
    y + p
  }
  y0 <- as.numeric(1:5)
  p <- 1
  i <- 0:10
  res <- difeq(y0, i, rhs, p, return_initial = TRUE)

  expect_equal(res, outer(y0, i, "+"))

  y0 <- runif(5)
  p <- runif(5)
  res <- difeq(y0, i, rhs, p, return_initial = TRUE)
  cmp <- y0 + outer(p, i)
  expect_equal(res, cmp)

  t0 <- runif(1)
  dt <- runif(1)
  res <- difeq(y0, i, rhs, p,
               return_initial = TRUE,
               n_history = length(i),
               t0 = t0, dt = dt)

  ## Check the history buffer:
  h <- attr(res, "history")
  expect_equal(h[1, ], i)
  expect_equal(h[2, ], t0 + i * dt)
  expect_equal(h[-(1:2), ], cmp)

  ## Check the returned values:
  attr(res, "history") <- NULL
  expect_equal(res, cmp)

  ## And again, but using the shortcut:
  res2 <- difeq(y0, max(i), rhs, p,
                return_initial = TRUE,
                n_history = length(i),
                t0 = t0, dt = dt)
  expect_equal(attr(res2, "history"), h)
  attr(res2, "history") <- NULL
  expect_equal(res2, res)

  ## And again, but with only a few history times:
  i2 <- seq(0, 10, by = 2)
  res <- difeq(y0, i2, rhs, p,
               return_initial = TRUE,
               n_history = length(i),
               t0 = t0, dt = dt)
  ## The history length should not change here.
  expect_identical(attr(res, "history"), h)
  attr(res, "history") <- NULL
  expect_equal(res, cmp[, i2 + 1])
})

test_that("output (R)", {
  rhs <- function(i, t, y, p) {
    ret <- y + p
    attr(ret, "output") <- sum(y)
    ret
  }
  y0 <- runif(5)
  p <- runif(5)
  i <- 0:10
  i2 <- seq(0, 10, by = 2)

  cmp <- y0 + outer(p, i)
  cmp2 <- y0 + outer(p, i2)

  res <- difeq(y0, i, rhs, p, return_initial = TRUE, n_out = 1)
  resh <- difeq(y0, i, rhs, p, return_initial = TRUE, n_out = 1,
                n_history = length(i))

  expect_equal(res, cmp, check.attributes = FALSE)
  expect_equal(resh, cmp, check.attributes = FALSE)

  output <- attr(res, "output")
  expect_is(res, "matrix")
  expect_equal(dim(output), c(1, length(i)))
  expect_equal(drop(output), colSums(res))
  expect_equal(attr(resh, "output"), output)

  res2 <- difeq(y0, i2, rhs, p, return_initial = TRUE, n_out = 1)
  resh2 <- difeq(y0, i2, rhs, p, return_initial = TRUE, n_out = 1,
                 n_history = length(i))

  expect_equal(res2, cmp2, check.attributes = FALSE)
  expect_equal(resh2, cmp2, check.attributes = FALSE)

  output2 <- attr(res2, "output")
  expect_is(res2, "matrix")
  expect_equal(dim(output2), c(1, length(i2)))
  expect_equal(drop(output2), colSums(res2))
  expect_equal(attr(resh2, "output"), output2)

  expect_equal(res, resh, check.attributes = FALSE)
  expect_equal(res2, resh2, check.attributes = FALSE)
})

test_that("logistic", {
  logistic <- function(i, t, y, p) {
    r * y * (1 - y)
  }
  y0 <- 0.1
  r <- 1.5
  i <- 0:20

  cmp <- difeq(y0, i, logistic, r, deSolve_compatible = TRUE)
  res <- difeq(y0, i, "logistic", r, dllname = "logistic",
               deSolve_compatible = TRUE)
  expect_equal(res, cmp)
})

test_that("error conditions", {
  rhs <- function(i, t, y, p) {
    y + p
  }
  y0 <- as.numeric(1:5)
  p <- 1
  i <- 0:10

  ## Beginning and end times are the same:
  expect_error(difeq(y0, 0, rhs, p),
               "Beginning and end times are the same")
  expect_error(difeq(y0, c(0, 0), rhs, p),
               "Beginning and end times are the same")
  expect_error(difeq(y0, c(i, 0), rhs, p),
               "Beginning and end times are the same")

  ## Incorrect order:
  expect_error(difeq(y0, rev(i), rhs, p),
               "Times not strictly increasing")
  expect_error(difeq(y0, c(0, 1, 1, 2), rhs, p),
               "Times not strictly increasing")
  expect_error(difeq(y0, c(0, 2, 1), rhs, p),
               "Times not strictly increasing")

  expect_error(difeq(y0, i, rhs, p, unknown = TRUE),
               "Invalid dot arguments")

  expect_error(difeq(y0, i, rhs, p, dllname = "logistic"),
               "dllname must not be given")

  expect_error(difeq(y0, (-5):(-1), rhs, p),
               "steps must be positive")
})
