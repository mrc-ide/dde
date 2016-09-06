check_ynames <- function(y, ynames, deSolve_compatible) {
  if (isTRUE(ynames)) {
    ynames <- names(y)
    ## This doesn't seem ideal but does produce more deSolve-like output
    if (deSolve_compatible && is.null(ynames)) {
      ynames <- as.character(seq_along(y))
    }
  } else if (is.null(ynames) || identical(as.vector(ynames), FALSE)) {
    ynames <- NULL
  } else if (is.character(ynames)) {
    if (length(ynames) != length(y)) {
      stop("ynames must be the same length as y")
    }
  } else {
    stop("Invalid value for ynames")
  }
  ynames
}

check_outnames <- function(n_out, outnames) {
  if (!is.null(outnames)) {
    if (is.character(outnames)) {
      if (length(outnames) != n_out) {
        stop("outnames must have length n_out")
      }
    } else {
      stop("Invalid value for outnames")
    }
  }
  outnames
}
