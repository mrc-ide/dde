context("shlib")

test_that("error reporting", {
  path <- tempfile(fileext = ".c")
  writeLines("not really c", path)
  expect_error(capture.output(shlib(path, "dde_")),
               "Error compiling source; see above for details")
})


test_that("don't rebuild dll", {
  expect_null(shlib("growth.c", "dde_"))
})
