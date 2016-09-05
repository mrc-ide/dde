context("interface")

test_that("invalid args", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  growth <- function(t, y, p) {
    y * p
  }

  expect_error(dopri(y, times, "lorenz", p, dllname = "lorenz",
                     unknown = 1),
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

  expect_error(dopri(y, times, growth, p, dllname = "lorenz",
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

  expect_error(dopri(y, times, "lorenz", p, dllname = "lorenz",
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

  expect_error(dopri(y, times, "lorenz", p, dllname = "lorenz",
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
  expect_identical(dopri5(y, times, "lorenz", p, dllname = "lorenz"),
                   dopri(y, times, "lorenz", p, dllname = "lorenz"))
  expect_identical(dopri5(y, times, "lorenz", p, dllname = "lorenz"),
                   dopri(y, times, "lorenz", p, dllname = "lorenz",
                         method = "dopri5"))
  expect_identical(dopri853(y, times, "lorenz", p, dllname = "lorenz"),
                   dopri(y, times, "lorenz", p, dllname = "lorenz",
                         method = "dopri853"))
})

test_that("invalid function input", {
  p <- c(10, 28, 8 / 3)
  y <- c(10, 1, 1)
  times <- seq(0, 1, length.out = 11)
  expect_error(dopri(y, times, 1, p),
               "Invalid input for 'fun'")
})

test_that("don't call ylag improperly", {
  expect_error(ylag(10), "Can't call this without being in an integration")
})
