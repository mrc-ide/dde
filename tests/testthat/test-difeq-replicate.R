context("discrete: replicated")

test_that("replication: simplest case", {
  rhs <- function(i, y, p) {
    y + p
  }

  y0 <- as.numeric(1:5)
  p <- 1
  i <- 0:10
  res <- difeq(y0, i, rhs, p)

  res2 <- difeq_replicate(3, y0, i, rhs, p)

  expect_is(res2, "array")
  expect_equal(dim(res2), c(11, 6, 3))
  expect_equal(res2[, , 1], res)
  expect_equal(res2[, , 2], res)
  expect_equal(res2[, , 3], res)
  expect_equal(names(attributes(res2)), "dim")
})


test_that("preserve y names", {
  rhs <- function(i, y, p) {
    y + p
  }

  y0 <- as.numeric(1:5)
  p <- 1
  i <- 0:10
  names(y0) <- letters[seq_along(y0)]
  res <- difeq(y0, i, rhs, p)

  res2 <- difeq_replicate(3, y0, i, rhs, p, as_array = TRUE)

  expect_equal(dimnames(res2), c(dimnames(res), list(NULL)))
})


test_that("varying initial conditions", {
  rhs <- function(i, y, p) {
    y + p
  }

  y0 <- lapply(1:3, seq, length.out = 5)
  p <- 1
  i <- 0:10
  cmp <- difeq(y0[[1]], i, rhs, p, return_step = FALSE)
  res <- difeq_replicate(3, y0, i, rhs, p, return_step = FALSE)

  expect_equal(res[, , 1L], cmp)
  expect_equal(res[, , 2L], cmp + 1)
  expect_equal(res[, , 3L], cmp + 2)
})


test_that("names from varying initial conditions", {
  rhs <- function(i, y, p) {
    y + p
  }

  y0 <- lapply(1:3, seq, length.out = 5)
  names(y0) <- c("X", "Y", "Z")
  for (i in seq_along(y0)) {
    names(y0[[i]]) <- letters[seq_len(5)]
  }
  p <- 1
  i <- 0:10
  cmp <- difeq(y0[[1]], i, rhs, p, return_step = FALSE)
  res3 <- difeq_replicate(3, y0, i, rhs, p,
                          return_step = FALSE, as_array = TRUE)

  expect_equal(dimnames(res3), c(dimnames(cmp), list(NULL)))
})


test_that("replicate invalid input", {
  rhs <- function(i, y, p) {
    y + p
  }
  y0 <- as.numeric(1:5)
  p <- 1
  i <- 0:10

  expect_error(difeq_replicate(0, y0, i, rhs, p, unknown = "whatever"),
               "Invalid dot arguments")
  expect_error(difeq_replicate(1, y0, i, rhs, p, as_array = FALSE),
               "as_array = FALSE no longer supported")
  expect_error(difeq_replicate(0, y0, i, rhs, p),
               "n must be positive")
  expect_error(difeq_replicate(2, list(y0), i, rhs, p),
               "'y' must have 2 elements")
  expect_error(difeq_replicate(2, list(y0, y0[-1]), i, rhs, p),
               "All 'y' lengths must be the same")
  expect_error(difeq_replicate(2, list(setNames(y0, letters[1:5]),
                                       setNames(y0, LETTERS[1:5])),
                               i, rhs, p),
               "All 'y' names must be the same")
  expect_error(difeq_replicate(1, y0, -1, rhs, p),
               "steps must be positive")
  expect_error(difeq_replicate(1, y0, i, rhs, p, n_history = 1),
               "If given, n_history must be at least 2")
  expect_error(difeq_replicate(1, y0, i, rhs, p, dllname = "whatever"),
               "dllname must not be given")
  expect_error(difeq_replicate(2, matrix(y0), i, rhs, p),
               "Expected '2' columns in 'y'")
})


test_that("basic replication", {
  rhs <- function(i, y, p) {
    y + p
  }

  y0 <- lapply(1:3, seq, length.out = 5)

  p <- 1
  i <- 0:10
  cmp <- difeq(y0[[1]], i, rhs, p)
  res <- difeq_replicate(3, y0, i, rhs, p)

  expect_equal(res[, , 1L], difeq(y0[[1L]], i, rhs, p))
  expect_equal(res[, , 2L], difeq(y0[[2L]], i, rhs, p))
  expect_equal(res[, , 3L], difeq(y0[[3L]], i, rhs, p))

  expect_identical(difeq_replicate(3, do.call(cbind, y0), i, rhs, p), res)
})


test_that("add output, time", {
  rhs <- function(i, y, p) {
    ret <- y + p
    attr(ret, "output") <- sum(y)
    ret
  }

  y0 <- list(
    as.numeric(1:5),
    as.numeric(2:6),
    as.numeric(3:7))
  p <- 3
  i <- 0:10

  ## This is what we're trying to achieve:
  a <- difeq(y0[[1]], i, rhs, p, n_out = 1L)
  b <- difeq(y0[[2]], i, rhs, p, n_out = 1L)
  c <- difeq(y0[[3]], i, rhs, p, n_out = 1L)
  d <- difeq_replicate(3, y0, i, rhs, p, n_out = 1L)

  expect_equal(dim(d), c(dim(a), 3))
  expect_equal(d[, , 1], a)
  expect_equal(d[, , 2], b)
  expect_equal(d[, , 3], c)

  e <- difeq_replicate(3, y0, i, rhs, p, n_out = 1L,
                       return_output_with_y = FALSE)
  expect_equal(attr(e, "output"), d[, 7, , drop = FALSE])

  attr(e, "output") <- NULL
  expect_equal(e, d[, 1:6, , drop = FALSE])
})


test_that("limitations of replication interface", {
  rhs <- function(i, y, p) {
    y + p
  }
  y0 <- as.numeric(1:5)
  p <- 1
  i <- 0:10
  expect_error(
    difeq_replicate(1, y0, i, rhs, p, restartable = TRUE),
    "Can't return pointer when n_replicates > 1")
  expect_error(
    difeq_replicate(1, y0, i, rhs, p, return_history = TRUE),
    "Can't return history when n_replicates > 1")
})


test_that("Minimal output", {
  rhs <- function(i, y, p) {
    ret <- y + p
    attr(ret, "output") <- sum(y)
    ret
  }

  y0 <- cbind(1:5, 2:6, 3:7)
  rownames(y0) <- letters[1:5]
  outnames <- "total"

  p <- 1
  cmp <- difeq(y0[[1]], 10, rhs, p)
  res <- difeq_replicate(3, y0, 10, rhs, p, n_out = 1L,
                         outnames = outnames,
                         return_minimal = TRUE)
  output <- attr(res, "output")

  expect_equal(dim(res), c(5, 10, 3))
  expect_equal(dim(output), c(1, 10, 3))
  expect_equal(rownames(res), rownames(y0))
  expect_equal(rownames(output), outnames)
})
