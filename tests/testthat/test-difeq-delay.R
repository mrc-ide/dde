context("discrete, with delays")

test_that("fib", {
  ## Here's a fun simple example; the Fibonacci sequence:
  ##
  ## TODO: Because I can't get the 0 into the history buffer to start
  ## with, this is not quite correct, so it's really fib(2) and up.
  fib <- function(i, y, p) {
    y + yprev(i - 1L)
  }
  y0 <- 1
  i <- 0:10
  res <- difeq(y0, i, fib, NULL, n_history = 2, return_step = FALSE)
  expect_equal(res[1:11], c(1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144))

  h <- attr(res, "history")
  expect_equal(attr(h, "n"), 1L)
  expect_equal(dim(h), c(2, 2))

  res2 <- difeq(y0, i, fib, NULL, return_step = FALSE, n_history = 20)
  h2 <- attr(res2, "history")
  expect_equal(dim(h2), c(2, 11))
  expect_equal(h2[, 10:11], h, check.attributes = FALSE)
})

test_that("prev and output", {
  growth <- function(i, y, p) {
    ret <- y + p
    attr(ret, "output") <- yprev(i - 1L)
    ret
  }

  n <- 5
  y0 <- runif(n)
  p <- runif(n)
  i <- 0:10
  i2 <- seq(0, 10, by = 2)

  cmp <- t(y0 + outer(p, i))
  cmp2 <- t(y0 + outer(p, i2))

  res <- difeq(y0, i, growth, p, return_step = FALSE, n_out = 5,
               n_history = 2L, return_output_with_y = FALSE)
  expect_equal(res, cmp, check.attributes = FALSE)
  output <- attr(res, "output")
  expect_equal(output, rbind(y0, cmp[-nrow(cmp), ], deparse.level = 0))

  res2 <- difeq(y0, i2, growth, p, return_step = FALSE, n_out = 5,
                n_history = 2L, return_output_with_y = FALSE)
  expect_equal(res2, cmp2, check.attributes = FALSE)
  output2 <- attr(res2, "output")
  expect_equal(output2, output[i2 + 1, ])
})

test_that("yprev permutations", {
  rhs <- function(i, y, p) {
    iprev <- i - 1
    if (p == "one") {
      for (i in seq_along(y)) {
        y[i] <- y[i] + yprev(iprev, i)
      }
    } else if (p == "idx") {
      y <- y + yprev(iprev, seq_along(y))
    } else { # all
      y <- y + yprev(iprev)
    }
    y
  }

  i <- seq(0:10)
  y0 <- runif(5)
  res1 <- difeq(y0, i, rhs, "one", return_initial = TRUE, n_history = 2L)
  res2 <- difeq(y0, i, rhs, "idx", return_initial = TRUE, n_history = 2L)
  res3 <- difeq(y0, i, rhs, "all", return_initial = TRUE, n_history = 2L)

  cmp <- matrix(y0, length(i), length(y0), byrow = TRUE)
  for (j in 2:length(i)) {
    cmp[j, ] <- cmp[j - 1, ] + cmp[max(1, j - 2), ]
  }

  expect_equal(res1[, -1], cmp, check.attributes = FALSE)
  expect_equal(res2, res1)
  expect_equal(res3, res1)
})

test_that("yprev invalid input", {
  rhs <- function(i, y, p) {
    type <- p$type
    lag <- p$lag
    idx <- p$idx
    iprev <- if (is.numeric(lag)) i - lag else lag
    if (type == "one") {
      for (i in seq_along(idx)) {
        j <- idx[[i]]
        yp <- yprev(iprev, j)
        y[j] <- y[j] + yp
      }
    } else if (type == "idx") {
      yp <- yprev(iprev, idx)
      y[idx] <- y[idx] + yp
    } else { # all
      y <- y + yprev(iprev)
    }
    y
  }

  i <- seq(0:10)
  y0 <- runif(5)
  cmp <- matrix(y0, length(y0), length(i))
  for (j in 2:length(i)) {
    cmp[, j] <- cmp[, j - 1] + cmp[, max(1, j - 2)]
  }
  cmp <- t(cmp)

  p <- list(lag = 1.0, idx = as.numeric(seq_along(y0)), type = "one")
  res1 <- difeq(y0, i, rhs, p, n_history = 2L, return_initial = TRUE)
  expect_equal(res1[, -1], cmp, check.attributes = FALSE)

  p <- list(lag = 1.0, idx = as.numeric(seq_along(y0)), type = "idx")
  res1 <- difeq(y0, i, rhs, p, n_history = 2L, return_initial = TRUE)
  expect_equal(res1[, -1], cmp, check.attributes = FALSE)

  p <- list(lag = 1.0, idx = as.numeric(seq_along(y0) + 1), type = "one")
  expect_error(difeq(y0, i, rhs, p, n_history = 2L, return_initial = TRUE),
               "out of bounds")
  p <- list(lag = 1.0, idx = as.integer(seq_along(y0) + 1L), type = "one")
  expect_error(difeq(y0, i, rhs, p, n_history = 2L, return_initial = TRUE),
               "out of bounds")

  p <- list(lag = 1.0, idx = "one", type = "one")
  expect_error(difeq(y0, i, rhs, p, n_history = 2L, return_initial = TRUE),
               "Invalid type")
  p <- list(lag = 1.0, idx = c("one", "two"), type = "idx")
  expect_error(difeq(y0, i, rhs, p, n_history = 2L, return_initial = TRUE),
               "Invalid type")

  p <- list(lag = - 1.0, type = "one", idx = 1L)
  expect_error(difeq(y0, i, rhs, p, n_history = 2L, return_initial = TRUE),
               "did not find step")
  p <- list(lag = - 1.0, type = "idx", idx = 1L)
  expect_error(difeq(y0, i, rhs, p, n_history = 2L, return_initial = TRUE),
               "did not find step")
  p <- list(lag = - 1.0, type = "all")
  expect_error(difeq(y0, i, rhs, p, n_history = 2L, return_initial = TRUE),
               "did not find step")

  p <- list(lag = NULL, type = "all")
  expect_error(difeq(y0, i, rhs, p, n_history = 2L, return_initial = TRUE),
               "Expected a scalar")
  p <- list(lag = c(1L, 2L), type = "all")
  expect_error(difeq(y0, i, rhs, p, n_history = 2L, return_initial = TRUE),
               "Expected a scalar")
  p <- list(lag = "one", type = "all")
  expect_error(difeq(y0, i, rhs, p, n_history = 2L, return_initial = TRUE),
               "Expected an integer")
})

test_that("integer lag", {
  growth <- function(i, y, p) {
    y + yprev(i - 1L)
  }

  i <- seq(0:10)
  y0 <- runif(5)
  p <- numeric(0) # no parameters

  cmp <- matrix(y0, length(y0), length(i))
  for (j in 2:length(i)) {
    cmp[, j] <- cmp[, j - 1] + cmp[, max(1, j - 2)]
  }
  cmp <- t(cmp)

  res <- difeq(y0, i, growth, p, n_history = 2L,
               return_history = FALSE, return_step = FALSE)
  expect_equal(res, cmp)

  res <- difeq(y0, i, "growth", p, n_history = 2L, dllname = "dde_growth_int",
               return_history = FALSE, return_step = FALSE)
  expect_equal(res, cmp)
})

test_that("restart", {
  growth <- function(i, y, p) {
    ret <- y + p
    attr(ret, "output") <- yprev(i - 5L)
    ret
  }

  n <- 5L
  y0 <- runif(n)
  p <- runif(n)

  tt <- 0:50
  tc <- 20
  tt1 <- tt[tt <= tc]
  tt2 <- tt[tt >= tc]

  cmp <- y0 + outer(p, tt)

  res <- difeq(y0, tt, growth, p, n_out = n, n_history = 100L,
               restartable = FALSE)
  h <- attr(res, "history")
  expect_null(attr(res, "ptr"))
  expect_null(attr(res, "restart_data"))

  res1 <- difeq(y0, tt1, growth, p,  n_out = n, n_history = 100L,
                restartable = TRUE)

  ## We have useful things here:
  expect_is(attr(res1, "ptr"), "externalptr")
  expect_is(attr(res1, "restart_data"), "list")

  ## Do the restart
  res2 <- difeq_continue(res1, tt2, copy = FALSE)
  h2 <- attr(res2, "history")

  ## This is key; do all history elements match up:
  expect_equal(h2, h)

  expect_equal(res2[1, ], res1[nrow(res1), ])
  expect_equal(res, rbind(res1, res2[-1, ]), check.attributes = FALSE)
})

test_that("change y on restart", {
  growth <- function(i, y, p) {
    ret <- y + p
    attr(ret, "output") <- yprev(i - 5L)
    ret
  }

  n <- 5L
  y0 <- runif(n)
  p <- runif(n)

  tt <- 0:50
  tc <- 20
  it2 <- tt >= tc
  tt1 <- tt[tt <= tc]
  tt2 <- tt[it2]

  i <- seq_len(n) + 1L
  j <- i + n

  res <- difeq(y0, tt, growth, p, n_out = n, n_history = 100L,
               restartable = FALSE)
  h <- attr(res, "history")

  res1 <- difeq(y0, tt1, growth, p,  n_out = n, n_history = 100L,
                restartable = TRUE)
  y1 <- res1[nrow(res1), i] + 1
  res2 <- difeq_continue(res1, tt2, y1, copy = FALSE)
  h2 <- attr(res2, "history")

  expect_equal(res2[, 1], res[it2, 1])
  ## All offset by one:
  expect_equal(res2[, i], res[it2, i] + 1)

  ## Need to work out what I did here to get the output confirmed;
  ## it's lagged, but can't recall by what.
  expect_equal(res[-(1:5), j],
               res[seq_len(nrow(res) - 5), i])
  ## Cool, we're off by one, which probably means that we failed to
  ## get the correct number out, or it's set incorrectly in the
  ## history (most likely).
  ##
  ## TODO: this would be easiest to debug with plain output rather
  ## than lagged.  More testing will be added and I can pick that up
  ## later.
  ##
  ##   expect_equal(res2[-(1:5), j],
  ##                res2[seq_len(nrow(res2) - 5), i])
})
