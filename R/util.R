assert_scalar <- function(x, name=deparse(substitute(x))) {
  if (length(x) != 1) {
    stop(sprintf("%s must be a scalar", name), call.=FALSE)
  }
}

assert_logical <- function(x, name=deparse(substitute(x))) {
  if (!is.logical(x)) {
    stop(sprintf("%s must be logical", name), call.=FALSE)
  }
}

assert_scalar_logical <- function(x, name=deparse(substitute(x))) {
  assert_scalar(x, name)
  assert_logical(x, name)
}

assert_size <- function(x, strict=FALSE, name=deparse(substitute(x))) {
  assert_scalar(x, name)
  assert_integer(x, strict, name)
  assert_nonnegative(x, name)
}

assert_nonnegative <- function(x, name=deparse(substitute(x))) {
  if (x < 0) {
    stop(sprintf("%s must be nonnegative", name), call.=FALSE)
  }
}

assert_integer <- function(x, strict=FALSE, name=deparse(substitute(x))) {
  if (!(is.integer(x))) {
    usable_as_integer <-
      !strict && is.numeric(x) && (max(abs(as.integer(x) - x)) < 1e-8)
    if (!usable_as_integer) {
      stop(sprintf("%s must be integer", name), call.=FALSE)
    }
  }
}
