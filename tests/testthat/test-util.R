context("util")

test_that("assertions work", {
  expect_error(assert_scalar(NULL), "must be a scalar")
  expect_error(assert_scalar(numeric(0)), "must be a scalar")
  expect_error(assert_scalar(1:2), "must be a scalar")

  expect_error(assert_scalar_logical(1), "must be logical")
  expect_error(assert_scalar_logical(NA), "must not be NA")
  expect_error(assert_scalar_logical(c(TRUE, TRUE)), "must be a scalar")

  expect_error(assert_size(1.5), "must be integer")
  expect_error(assert_size(-2L), "must be nonnegative")
  expect_error(assert_size(NA_integer_), "must not be NA")
  expect_error(assert_size(c(1, 2)), "must be a scalar")

  expect_error(assert_scalar_character(character(0)), "must be a scalar")
  expect_error(assert_scalar_character(c("a", "b")), "must be a scalar")
  expect_error(assert_scalar_character(1), "must be character")
  expect_error(assert_scalar_character(NA_character_), "must not be NA")

  expect_error(match_value("a", c("b", "c")),
               "must be one of {b, c}", fixed = TRUE)

  expect_error(assert_numeric("a"), "must be numeric")
  expect_error(assert_numeric(TRUE), "must be numeric")
  expect_error(assert_scalar_numeric(NA_real_), "must not be NA")
  expect_error(assert_scalar_numeric(c(1, 2)), "must be a scalar")
  expect_silent(assert_numeric(1L))

  expect_error(assert_positive(0L), "must be positive")
  expect_silent(assert_positive(1L))
})
