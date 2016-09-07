context("discrete, with delays")

test_that("fib", {
  ## Here's a fun simple example; the Fibonacci sequence:
  ##
  ## TODO: Because I can't get the 0 into the history buffer to start
  ## with, this is not quite correct, so it's really fib(2) and up.
  fib <- function(i, t, y, p) {
    y + yprev(i - 1L)
  }
  y0 <- 1
  i <- 0:10
  res <- difeq(y0, i, fib, NULL, return_initial = TRUE, n_history = 2)
  expect_equal(res[1:11], c(1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144))

  h <- attr(res, "history")
  expect_equal(attr(h, "n"), 1L)
  expect_equal(dim(h), c(3, 2))

  res2 <- difeq(y0, i, fib, NULL, return_initial = TRUE, n_history = 20)
  h2 <- attr(res2, "history")
  expect_equal(dim(h2), c(3, 11))
  expect_equal(h2[, 10:11], h, check.attributes = FALSE)
})

test_that("prev and output", {
  growth <- function(i, t, y, p) {
    ret <- y + p
    attr(ret, "output") <- yprev(i - 1L)
    ret
  }

  y0 <- runif(5)
  p <- runif(5)
  i <- 0:10
  i2 <- seq(0, 10, by = 2)

  cmp <- y0 + outer(p, i)
  cmp2 <- y0 + outer(p, i2)

  res <- difeq(y0, i, growth, p, return_initial = TRUE, n_out = 5,
               n_history = 2L)
  expect_equal(res, cmp, check.attributes = FALSE)
  output <- attr(res, "output")
  expect_equal(output, cbind(y0, cmp[, -ncol(cmp)], deparse.level = 0))

  res2 <- difeq(y0, i2, growth, p, return_initial = TRUE, n_out = 5,
                n_history = 2L)
  expect_equal(res2, cmp2, check.attributes = FALSE)
  output2 <- attr(res2, "output")
  expect_equal(output2, output[, i2 + 1])
})
