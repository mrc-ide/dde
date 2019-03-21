context("discrete")

## Solve some basic discrete time equations
test_that("increase", {
  rhs <- function(i, y, p) {
    y + p
  }
  y0 <- as.numeric(1:5)
  p <- 1
  i <- 0:10
  res <- difeq(y0, i, rhs, p)

  expect_equal(res[, 1], i)
  expect_equal(unname(res[, -1]), outer(i, y0, "+"))

  y0 <- runif(5)
  p <- runif(5)
  res <- difeq(y0, i, rhs, p)
  cmp <- t(y0 + outer(p, i))
  expect_equal(unname(res[, -1]), cmp)

  res <- difeq(y0, i, rhs, p, n_history = length(i))

  ## Check the history buffer:
  h <- attr(res, "history")
  expect_equal(h[1, ], i)
  expect_equal(h[-1, ], t(cmp))

  ## Check the returned values:
  attr(res, "history") <- NULL
  expect_equal(unname(res[, -1]), cmp)

  ## And again, but using the shortcut:
  res2 <- difeq(y0, max(i), rhs, p,
                return_initial = TRUE,
                n_history = length(i))
  expect_equal(attr(res2, "history"), h)
  attr(res2, "history") <- NULL
  expect_equal(res2, res)

  ## And again, but with only a few history times:
  i2 <- seq(0, 10, by = 2)
  res <- difeq(y0, i2, rhs, p,
               return_initial = TRUE,
               n_history = length(i))
  ## The history length should not change here.
  expect_identical(attr(res, "history"), h)
  attr(res, "history") <- NULL
  expect_equal(unname(res[, -1]), cmp[i2 + 1, ])
})

test_that("output (R)", {
  rhs <- function(i, y, p) {
    ret <- y + p
    attr(ret, "output") <- sum(y)
    ret
  }
  y0 <- runif(5)
  p <- runif(5)
  i <- 0:10
  i2 <- seq(0, 10, by = 2)

  cmp <- y0 + outer(p, i)
  cmp2 <- cmp[, i2 + 1L]

  ## We run this three times:
  ##   * res_o: coming off the output storage
  ##   * res_h: coming off history
  ##   * res_n: coming off of short history
  res_o <- difeq(y0, i, rhs, p, return_by_column = FALSE, return_step = FALSE,
                 return_output_with_y = FALSE,
                 n_out = 1)
  res_h <- difeq(y0, i, rhs, p, return_by_column = FALSE, return_step = FALSE,
                 return_output_with_y = FALSE,
                 n_out = 1, n_history = length(i))
  res_n <- difeq(y0, i, rhs, p, return_by_column = FALSE, return_step = FALSE,
                 return_output_with_y = FALSE,
                 n_out = 1, n_history = 2L)

  expect_equal(res_o, cmp, check.attributes = FALSE)
  expect_equal(res_h, cmp, check.attributes = FALSE)
  expect_equal(res_n, cmp, check.attributes = FALSE)

  output <- attr(res_o, "output")
  expect_is(output, "matrix")
  expect_equal(dim(output), c(1, length(i)))
  ## This checks off-by-one errors
  expect_equal(drop(output), colSums(res_o))
  expect_equal(drop(attr(res_h, "output")), colSums(res_h))
  expect_equal(drop(attr(res_n, "output")), colSums(res_n))

  ## Also check with no initial condition; this changes things around
  ## a bit.
  res_o <- difeq(y0, i, rhs, p, return_initial = FALSE,
                 return_by_column = FALSE, return_step = FALSE,
                 return_output_with_y = FALSE, n_out = 1)
  res_h <- difeq(y0, i, rhs, p, return_initial = FALSE,
                 return_by_column = FALSE, return_step = FALSE,
                 return_output_with_y = FALSE, n_out = 1,
                 n_history = length(i))
  res_n <- difeq(y0, i, rhs, p, return_initial = FALSE,
                 return_by_column = FALSE, return_step = FALSE,
                 return_output_with_y = FALSE, n_out = 1,
                 n_history = 2L)

  expect_equal(res_o, cmp[, -1], check.attributes = FALSE)
  expect_equal(res_h, cmp[, -1], check.attributes = FALSE)
  expect_equal(res_n, cmp[, -1], check.attributes = FALSE)

  expect_equal(drop(attr(res_o, "output")), colSums(res_o))
  expect_equal(drop(attr(res_h, "output")), colSums(res_h))
  expect_equal(drop(attr(res_n, "output")), colSums(res_n))

  ## There's another huge class of bugs that turns up when output is
  ## filtered, so try that too.
  res_o <- difeq(y0, i2, rhs, p, return_by_column = FALSE, return_step = FALSE,
                 return_output_with_y = FALSE, n_out = 1)
  res_h <- difeq(y0, i2, rhs, p, return_by_column = FALSE, return_step = FALSE,
                 return_output_with_y = FALSE, n_out = 1,
                 n_history = length(i))
  res_n <- difeq(y0, i2, rhs, p, return_by_column = FALSE, return_step = FALSE,
                 return_output_with_y = FALSE, n_out = 1,
                 n_history = 2L)

  expect_equal(res_o, cmp2, check.attributes = FALSE)
  expect_equal(res_h, cmp2, check.attributes = FALSE)
  expect_equal(res_n, cmp2, check.attributes = FALSE)

  output <- attr(res_o, "output")
  expect_is(output, "matrix")
  expect_equal(dim(output), c(1, length(i2)))
  expect_equal(drop(output), colSums(res_o))
  expect_equal(drop(output), colSums(res_o))
  expect_equal(drop(attr(res_h, "output")), colSums(res_h))
  expect_equal(drop(attr(res_n, "output")), colSums(res_n))

  ## And while dropping initial conditions:
  res_o <- difeq(y0, i2, rhs, p, return_initial = FALSE,
                 return_by_column = FALSE, return_step = FALSE,
                 return_output_with_y = FALSE, n_out = 1)
  res_h <- difeq(y0, i2, rhs, p, return_initial = FALSE,
                 return_by_column = FALSE, return_step = FALSE,
                 return_output_with_y = FALSE, n_out = 1,
                 n_history = length(i))
  res_n <- difeq(y0, i2, rhs, p, return_initial = FALSE,
                 return_by_column = FALSE, return_step = FALSE,
                 return_output_with_y = FALSE, n_out = 1,
                 n_history = 2L)

  expect_equal(res_o, cmp2[, -1], check.attributes = FALSE)
  expect_equal(res_h, cmp2[, -1], check.attributes = FALSE)
  expect_equal(res_n, cmp2[, -1], check.attributes = FALSE)

  expect_equal(drop(attr(res_o, "output")), colSums(res_h))
  expect_equal(drop(attr(res_h, "output")), colSums(res_h))
  expect_equal(drop(attr(res_n, "output")), colSums(res_n))
})

test_that("transpose output", {
  rhs <- function(i, y, p) {
    ret <- y + p
    attr(ret, "output") <- sum(y)
    ret
  }
  y0 <- runif(5)
  p <- runif(5)
  i <- 0:10
  i2 <- seq(0, 10, by = 2)

  cmp <- y0 + outer(p, i)
  cmp2 <- cmp[, i2 + 1L]

  res_mo <- difeq(y0, i, rhs, p, return_minimal = TRUE, n_out = 1)
  res_mh <- difeq(y0, i, rhs, p, return_minimal = TRUE, n_out = 1,
                 n_history = length(i), return_history = FALSE)
  res_fo <- difeq(y0, i, rhs, p, n_out = 1, ynames = TRUE)
  res_fh <- difeq(y0, i, rhs, p, n_out = 1, ynames = TRUE,
                 n_history = length(i), return_history = FALSE)

  expect_equal(colnames(res_fo), c("step", as.character(seq_along(y0)), ""))
  expect_null(dimnames(res_mo))

  expect_equal(res_fo[, 1], i)
  expect_equal(dim(res_fo), c(ncol(res_mo) + 1, nrow(res_mo) + 2))

  expect_equal(res_mo, cmp[, -1], check.attributes = FALSE)
  expect_equal(res_fo[, -c(1, ncol(res_fo))], t(cmp), check.attributes = FALSE)

  expect_equal(res_fo[, ncol(res_fo)], colSums(cmp))
  expect_equal(attr(res_mo, "output")[1, ], colSums(cmp)[-1])

  expect_null(attr(res_fo, "output"))
  expect_equal(res_fo, res_fh)
  expect_equal(res_mo, res_mh)
})

test_that("incorrect output", {
  f <- function(i, y, p) {
    structure(y, output = p)
  }
  y0 <- runif(5)
  i <- 0:10

  ## Check that things generally work:
  expect_equal(difeq(y0, i, f, NULL, return_minimal = TRUE),
               matrix(y0, length(y0), length(i) - 1L))

  ## These are going to cause leaks for now, but I need to handle
  ## errors in the calling function too.
  ##
  ## The options for doing this are going to be either wrapping every
  ## call to an R objective function in tryCatch (which is quite slow;
  ## being ~3us rather than a few seconds) or smart pointers.
  expect_error(difeq(y0, i, f, NULL, n_out = 1),
               "Missing output")
  expect_error(difeq(y0, i, f, c(1, 2), n_out = 1),
               "Incorrect length")
  expect_error(difeq(y0, i, f, 1L, n_out = 1),
               "Incorrect type")
})

test_that("logistic", {
  logistic <- function(i, y, p) {
    r * y * (1 - y)
  }
  y0 <- 0.1
  r <- 1.5
  i <- 0:20

  cmp <- difeq(y0, i, logistic, r)
  res <- difeq(y0, i, "logistic", r, dllname = "dde_logistic")
  expect_equal(res, cmp)
})

test_that("vector output (R)", {
  growth <- function(i, y, p) {
    structure(y + p, output = y + 1)
  }

  y0 <- runif(5)
  p <- runif(5)
  i <- 0:10
  i2 <- seq(0, 10, by = 2)

  cmp <- y0 + outer(p, i)
  cmp2 <- y0 + outer(p, i2)

  ## TODO: The NA in history here looks very fixable to me.
  res <- difeq(y0, i, growth, p, return_initial = TRUE, n_out = 5L,
               return_by_column = FALSE, return_step = FALSE,
               return_output_with_y = FALSE,
               n_history = length(i))
  expect_equal(res, cmp, check.attributes = FALSE)
  expect_equal(attr(res, "output"), cmp + 1)
  h <- attr(res, "history")
  expect_equal(h[1, ], i)
  expect_equal(h[seq_along(y0) + 1, ], cmp)
  expect_equal(h[seq_along(y0) + 1 + length(y0), ],
               cbind(NA, cmp[, -1], deparse.level = 0) + 1)

  res2 <- difeq(y0, i2, growth, p, return_initial = TRUE, n_out = 5L,
               return_by_column = FALSE, return_step = FALSE,
               return_output_with_y = FALSE,
                n_history = length(i))
  j <- seq_len(length(y0) + 1)
  expect_equal(attr(res2, "history")[j, ], h[j, ])

  expect_equal(res2, cmp2, check.attributes = FALSE)
  expect_equal(attr(res2, "output"), cmp2 + 1)

  ## NOTE: At the moment the status of output within history is a bit
  ## of a mess; nobody should rely on it being there, or use it for
  ## anything.
})

test_that("error conditions", {
  rhs <- function(i, y, p) {
    y + p
  }
  y0 <- as.numeric(1:5)
  p <- 1
  i <- 0:10

  ## Beginning and end times are the same:
  expect_error(difeq(y0, 0, rhs, p),
               "Beginning and end steps are the same")
  expect_error(difeq(y0, c(0, 0), rhs, p),
               "Beginning and end steps are the same")
  expect_error(difeq(y0, c(i, 0), rhs, p),
               "Beginning and end steps are the same")

  ## Incorrect order:
  expect_error(difeq(y0, rev(i), rhs, p),
               "Steps not strictly increasing")
  expect_error(difeq(y0, c(0, 1, 1, 2), rhs, p),
               "Steps not strictly increasing")
  expect_error(difeq(y0, c(0, 2, 1), rhs, p),
               "Steps not strictly increasing")

  expect_error(difeq(y0, i, rhs, p, unknown = TRUE),
               "Invalid dot arguments")

  res <- difeq(y0, i, rhs, p, restartable = TRUE)
  expect_error(difeq_continue(res, i + 1, unknown = TRUE),
               "Invalid dot arguments")

  expect_error(difeq(y0, i, rhs, p, dllname = "dde_logistic"),
               "dllname must not be given")

  expect_error(difeq(y0, (-5):(-1), rhs, p),
               "steps must be positive")

  expect_error(difeq(y0, i, rhs, p, n_history = 1L),
               "n_history must be at least 2")
  expect_error(difeq(y0, i, rhs, p, n_history = -1L),
               "n_history must be nonnegative")
})

test_that("names", {
  p <- c(10, 28, 8 / 3)
  y0 <- c(10, 1, 1)

  rhs <- function(i, y, p) {
    structure(y + p$x, output = p$output)
  }

  y0 <- as.numeric(1:5)
  p <- list(x = 1, output = NULL)
  i <- 0:10

  nms <- letters[seq_along(y0)]
  cmp <- list(NULL, c("step", nms))

  ## These are duplicated from ode:
  expect_null(dimnames(difeq(y0, i, rhs, p, ynames = FALSE)))
  expect_equal(dimnames(difeq(y0, i, rhs, p, ynames = nms)), cmp)
  expect_equal(dimnames(difeq(setNames(y0, nms), i, rhs, p)), cmp)
  expect_null(dimnames(difeq(setNames(y0, nms), i, rhs, p, ynames = FALSE)),
              cmp)

  p$output <- runif(3)
  onms <- LETTERS[seq_along(p$output)]
  ocmp <- list(NULL, onms)

  f <- function(..., return_output_with_y = FALSE) {
    difeq(y0, i, rhs, p, n_out = length(p$output),
          return_output_with_y = return_output_with_y, ...)
  }
  expect_equal(dim(attr(f(), "output")), c(length(i), length(p$output)))
  expect_null(dimnames(attr(f(), "output")))
  expect_null(dimnames(attr(f(outnames = NULL), "output")))
  expect_equal(dimnames(attr(f(outnames = onms), "output")), ocmp)
  expect_error(f(outnames = nms), "outnames must have length n_out")
  expect_error(f(outnames = 1), "Invalid value for outnames")

  ## Check both together:
  res <- f(ynames = nms, outnames = onms)
  expect_equal(dimnames(res), cmp)
  expect_equal(dimnames(attr(res, "output")), ocmp)
})

test_that("don't call yprev improperly", {
  expect_error(yprev(10), "Can't call this without being in an integration")
})

test_that("externalptr input", {
  y0 <- runif(5)
  r <- runif(length(y0))
  i <- 0:11

  cmp <- difeq(y0, i, "logistic", r, dllname = "dde_logistic")

  ptr <- .Call("logistic_init",  r, PACKAGE = "dde_logistic2")
  expect_is(ptr, "externalptr")
  res <- difeq(y0, i, "logistic", ptr, parms_are_real = FALSE,
               dllname = "dde_logistic2")

  expect_identical(res, cmp)
})

test_that("externalptr input safety", {
  y0 <- runif(5)
  r <- runif(length(y0))
  i <- 0:11
  ptr <- .Call("logistic_init",  r, PACKAGE = "dde_logistic2")
  expect_error(difeq(y0, i, "logistic", make_null_pointer(ptr),
                     parms_are_real = FALSE,
                     dllname = "dde_logistic2"),
               "Was passed null pointer for 'parms'")
})

test_that("externalptr target", {
  y0 <- runif(5)
  r <- runif(length(y0))
  i <- 0:11
  cmp <- difeq(y0, i, "logistic", r, dllname = "dde_logistic")
  ptr <- getNativeSymbolInfo("logistic", PACKAGE = "dde_logistic")
  expect_identical(difeq(y0, i, ptr, r), cmp)
})

test_that("externalptr target safety", {
  y0 <- runif(5)
  r <- runif(length(y0))
  i <- 0:11
  ptr <- getNativeSymbolInfo("logistic", PACKAGE = "dde_logistic")
  expect_error(difeq(y0, i, make_null_pointer(ptr), r),
               "Was passed null pointer for 'target'")
})

test_that("grow history", {
  rhs <- function(i, y, p) {
    y + p
  }
  y0 <- runif(5)
  p <- 1
  i <- 0:50

  res <- difeq(y0, i, rhs, p, return_initial = FALSE,
               n_history = 5, grow_history = TRUE)
  h <- attr(res, "history")
  expect_equal(ncol(h), length(i))

  cmp <- difeq(y0, i, rhs, p, return_initial = FALSE, n_history = 100)
  hc <- attr(cmp, "history")

  expect_equal(h, hc)
  expect_equal(cmp, res)
})


test_that("restart and copy", {
  growth <- function(i, y, p) {
    y + p
  }

  n <- 5L
  y0 <- runif(n)
  p <- runif(n)

  tt <- 0:50
  tc <- 20
  tt1 <- tt[tt <= tc]
  tt2 <- tt[tt >= tc]

  cmp <- difeq(y0, tt, growth, p)
  cmp2 <- cmp[tt >= tc, ]

  res1 <- difeq(y0, tt1, growth, p, restartable = TRUE)

  res2 <- difeq_continue(res1, tt2, copy = TRUE)
  res3 <- difeq_continue(res1, tt2, copy = TRUE)

  expect_equal(res2, cmp2, check.attributes = FALSE)
  expect_equal(res2, res3)

  res4 <- difeq_continue(res1, tt2, copy = FALSE)
  expect_equal(res4, cmp2, check.attributes = FALSE)

  ## Can't step forward now:
  expect_error(difeq_continue(res1, tt2),
               "Incorrect initial step on simulation restart")
})

test_that("restart error handling", {
  growth <- function(i, y, p) {
    y + p
  }
  n <- 5L
  y0 <- runif(n)
  p <- runif(n)
  tt1 <- 0:10
  tt2 <- 10:20
  res1 <- difeq(y0, tt1, growth, p, restartable = TRUE)

  y1 <- res1[nrow(res1), seq_len(n) + 1]
  expect_error(difeq_continue(res1, tt2, y1[-1], copy = FALSE),
               "Incorrect size 'y' on simulation restart", fixed = TRUE)
  expect_error(difeq_continue(res1, tt2[1], y1, copy = FALSE),
               "At least two steps must be given", fixed = TRUE)
  expect_error(difeq_continue(res1, numeric(0), y1, copy = FALSE),
               "At least two steps must be given", fixed = TRUE)
  expect_error(difeq_continue(res1, tt2[-1], y1, copy = FALSE),
               "Incorrect initial step on simulation restart", fixed = TRUE)
})

test_that("incorrect length output", {
  f <- function(i, y, p) {
    1
  }
  expect_error(difeq(c(1, 2), 0:5, f, NULL),
               "Incorrect length variable output")
})
