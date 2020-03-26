context("discrete: batch")

test_that("replication: simplest case", {
  rhs <- function(i, y, p) {
    y + p
  }

  y0 <- as.numeric(1:5)
  p <- as.list(1:3)
  i <- 0:10

  ## This is what we're trying to achieve:
  a <- difeq(y0, i, rhs, p[[1]], return_step = FALSE)
  b <- difeq(y0, i, rhs, p[[2]], return_step = FALSE)
  c <- difeq(y0, i, rhs, p[[3]], return_step = FALSE)
  d <- difeq_replicate(3, y0, i, rhs, p[[1]], return_step = FALSE)

  res <- difeq_batch(y0, i, rhs, p, return_step = FALSE)

  expect_equal(dim(res), dim(d))
  expect_equal(res[, , 1], a)
  expect_equal(res[, , 2], b)
  expect_equal(res[, , 3], c)
})


test_that("add output, time", {
  rhs <- function(i, y, p) {
    ret <- y + p
    attr(ret, "output") <- sum(y)
    ret
  }

  y0 <- as.numeric(1:5)
  p <- as.list(1:3)
  i <- 0:10

  ## This is what we're trying to achieve:
  a <- difeq(y0, i, rhs, p[[1]], n_out = 1L)
  b <- difeq(y0, i, rhs, p[[2]], n_out = 1L)
  c <- difeq(y0, i, rhs, p[[3]], n_out = 1L)
  d <- difeq_replicate(3, y0, i, rhs, p[[1]], n_out = 1L)

  res <- difeq_batch(y0, i, rhs, p, n_out = 1L)

  expect_equal(dim(res), dim(d))
  expect_equal(res[, , 1], a)
  expect_equal(res[, , 2], b)
  expect_equal(res[, , 3], c)
})
