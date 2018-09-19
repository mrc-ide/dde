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


test_that("replication: don't simplify", {
  rhs <- function(i, y, p) {
    y + p
  }

  y0 <- as.numeric(1:5)
  p <- 1
  i <- 0:10
  res <- difeq(y0, i, rhs, p)

  res2 <- difeq_replicate(3, y0, i, rhs, p, as_array = FALSE)
  expect_equal(res2, rep(list(res), 3))
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
  res3 <- difeq_replicate(3, y0, i, rhs, p, as_array = FALSE)

  expect_equal(dimnames(res2), c(dimnames(res), list(NULL)))
  expect_equal(res3, rep(list(res), 3))
})


test_that("restartable", {
  rhs <- function(i, y, p) {
    y + p
  }

  y0 <- as.numeric(1:5)
  p <- 1
  i <- 0:10
  names(y0) <- letters[seq_along(y0)]
  res <- difeq(y0, i, rhs, p, restartable = TRUE)

  res2 <- difeq_replicate(3, y0, i, rhs, p,
                          restartable = TRUE, as_array = FALSE)
  res3 <- difeq_replicate(3, y0, i, rhs, p,
                          restartable = TRUE, as_array = TRUE)

  expect_is(attr(res2[[1]], "restart_data"), "list")
  expect_is(attr(res3, "restart_data"), "list")
  expect_equal(names(attr(res2[[1]], "restart_data")),
               names(attr(res, "restart_data")))
  expect_equal(names(attr(res3, "restart_data")[[1]]),
               names(attr(res, "restart_data")))

  expect_is(attr(res2[[1]], "ptr"), "externalptr")
  expect_is(attr(res3, "ptr")[[1]], "externalptr")
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
  res2 <- difeq_replicate(3, y0, i, rhs, p,
                          return_step = FALSE, as_array = FALSE)
  res3 <- difeq_replicate(3, y0, i, rhs, p,
                          return_step = FALSE, as_array = TRUE)

  expect_equal(names(res2), names(y0))
  expect_equal(dimnames(res2[[1]]), dimnames(cmp))
  expect_equal(dimnames(res2[[3]]), dimnames(cmp))

  expect_equal(dimnames(res3), c(dimnames(cmp), list(names(y0))))
})


test_that("replicate invalid input", {
  rhs <- function(i, y, p) {
    y + p
  }
  y0 <- as.numeric(1:5)
  p <- 1
  i <- 0:10

  expect_error(difeq_replicate(0, y0, i, rhs, p),
               "n must be positive")
  expect_error(difeq_replicate(2, list(y0), i, rhs, p),
               "'y' must have 2 elements")
  expect_error(difeq_replicate(2, list(y0, y0[-1]), i, rhs, p),
               "All 'y' lengths must be the same")
})
