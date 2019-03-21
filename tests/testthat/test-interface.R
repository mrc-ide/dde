context("interface")

test_that("invalid args", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  growth <- function(t, y, p) {
    y * p
  }

  expect_error(dopri(y, times, "lorenz", p, dllname = "dde_lorenz",
                     unknown = 1),
               "Invalid dot arguments")

  res <- dopri(y, times, "lorenz", p, dllname = "dde_lorenz",
               restartable = TRUE)
  times2 <- seq(1, 2, length.out = 11)
  expect_error(dopri_continue(res, times2, unknown = 1),
               "Invalid dot arguments")
})

test_that("R function with dllname", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  growth <- function(t, y, p) {
    y * p
  }
  total <- function(t, y, p) {
    sum(y)
  }

  expect_error(dopri(y, times, growth, p, dllname = "dde_lorenz",
                     n_out = 1, output = total),
               "dllname must not be given")
})

test_that("R function with dllname", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  growth <- function(t, y, p) {
    y * p
  }
  total <- function(t, y, p) {
    sum(y)
  }

  expect_error(dopri(y, times, "lorenz", p, dllname = "dde_lorenz",
                     n_out = 1, output = total),
               "output must be a compiled function")
  expect_error(dopri(y, times, growth, p,
                     n_out = 2L, output = "lorenz_output"),
               "output must be an R function")
})

test_that("output with no n_out", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  growth <- function(t, y, p) {
    y * p
  }
  total <- function(t, y, p) {
    sum(y)
  }

  expect_error(dopri(y, times, "lorenz", p, dllname = "dde_lorenz",
                     output = "lorenz_output"),
               "n_out must be specified")
  expect_error(dopri(y, times, growth, p,
                     output = total),
               "n_out must be specified")
})

test_that("helpers", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  expect_identical(dopri5(y, times, "lorenz", p, dllname = "dde_lorenz"),
                   dopri(y, times, "lorenz", p, dllname = "dde_lorenz"))
  expect_identical(dopri5(y, times, "lorenz", p, dllname = "dde_lorenz"),
                   dopri(y, times, "lorenz", p, dllname = "dde_lorenz",
                         method = "dopri5"))
  expect_identical(dopri853(y, times, "lorenz", p, dllname = "dde_lorenz"),
                   dopri(y, times, "lorenz", p, dllname = "dde_lorenz",
                         method = "dopri853"))
})

test_that("invalid function input", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  expect_error(dopri(y, times, 1, p),
               "Invalid input for 'func'")
})

test_that("don't call ylag improperly", {
  expect_error(ylag(10), "Can't call this without being in an integration")
})

test_that("Missing output function", {
  growth <- function(t, y, p) y
  output <- function(t, y, p) sum(y)
  expect_error(dopri(1, 0:10, growth, NULL, n_out = 1),
               "Invalid input for 'output'")
})

test_that("Zero history in a lag model", {
  growth <- function(t, y, p) y
  output <- function(t, y, p) ylag(t - 2.0)
  tt <- seq(0, 10, length.out = 101)
  expect_error(dopri(1, tt, growth, NULL, n_out = 1, output = output),
               "Integration failure: can't use ylag in model with no history")
})
