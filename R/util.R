assert_scalar <- function(x, name = deparse(substitute(x))) {
  if (length(x) != 1) {
    stop(sprintf("%s must be a scalar", name), call. = FALSE)
  }
}


assert_logical <- function(x, name = deparse(substitute(x))) {
  if (!is.logical(x)) {
    stop(sprintf("%s must be logical", name), call. = FALSE)
  }
}


assert_nonmissing <- function(x, name = deparse(substitute(x))) {
  if (any(is.na(x))) {
    stop(sprintf("%s must not be NA", name), call. = FALSE)
  }
}


assert_scalar_logical <- function(x, name = deparse(substitute(x))) {
  assert_scalar(x, name)
  assert_logical(x, name)
  assert_nonmissing(x, name)
}


assert_character <- function(x, name = deparse(substitute(x))) {
  if (!is.character(x)) {
    stop(sprintf("%s must be character", name), call. = FALSE)
  }
}


assert_scalar_character <- function(x, name = deparse(substitute(x))) {
  assert_scalar(x, name)
  assert_character(x, name)
  assert_nonmissing(x, name)
}


assert_numeric <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x)) {
    stop(sprintf("%s must be numeric", name), call. = FALSE)
  }
}


assert_scalar_numeric <- function(x, name = deparse(substitute(x))) {
  assert_scalar(x, name)
  assert_numeric(x, name)
  assert_nonmissing(x, name)
}


assert_size <- function(x, strict = FALSE, name = deparse(substitute(x))) {
  assert_scalar(x, name)
  assert_integer(x, strict, name)
  assert_nonmissing(x, name)
  assert_nonnegative(x, name)
}


assert_positive_integer <- function(x, strict = FALSE,
                                    name = deparse(substitute(x))) {
  assert_scalar(x, name)
  assert_integer(x, strict, name)
  assert_nonmissing(x, name)
  assert_positive(x, name)
}


assert_nonnegative <- function(x, name = deparse(substitute(x))) {
  if (x < 0) {
    stop(sprintf("%s must be nonnegative", name), call. = FALSE)
  }
}


assert_positive <- function(x, name = deparse(substitute(x))) {
  if (x <= 0) {
    stop(sprintf("%s must be positive", name), call. = FALSE)
  }
}


assert_integer <- function(x, strict = FALSE, name = deparse(substitute(x))) {
  if (!(is.integer(x))) {
    usable_as_integer <-
      !strict && is.numeric(x) && (max(abs(as.integer(x) - x)) < 1e-8)
    if (!usable_as_integer) {
      stop(sprintf("%s must be integer", name), call. = FALSE)
    }
  }
}


assert_length <- function(x, n, name = deparse(substitute(x))) {
  if (length(x) != n) {
    stop(sprintf("'%s' must have %d elements (recieved %d)",
                 name, n, length(x)), call. = FALSE)
  }
}


match_value <- function(x, choices, name = deparse(substitute(x))) {
  assert_scalar_character(x, name)
  i <- match(x, choices)
  if (is.na(i)) {
    stop(sprintf("%s must be one of {%s}", name,
                 paste(choices, collapse = ", ")),
         call. = FALSE)
  }
  choices[[i]]
}


logopt <- function(value, default, name = deparse(substitute(value))) {
  if (is.null(value)) {
    default
  } else {
    assert_scalar_logical(value, name)
    value
  }
}


is_native_symbol_info <- function(x) {
  inherits(x, "NativeSymbolInfo")
}


collect_attributes <- function(object, name, target) {
  if (!is.null(attr(object[[1L]], name, exact = TRUE))) {
    attr(target, name) <- lapply(object, "attr", name, exact = TRUE)
  }
  target
}
