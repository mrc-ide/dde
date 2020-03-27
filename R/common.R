check_ynames <- function(y, ynames) {
  if (isTRUE(ynames)) {
    ynames <- if (is.matrix(y)) rownames(y) else names(y)
    ## This doesn't seem ideal but does produce more deSolve-like output
    if (is.null(ynames)) {
      ynames <- as.character(seq_along(y))
    }
  } else if (is.null(ynames) || identical(as.vector(ynames), FALSE)) {
    ynames <- NULL
  } else if (is.character(ynames)) {
    len <- if (is.matrix(y)) nrow(y) else length(y)
    if (length(ynames) != len) {
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


prepare_output <- function(ret, times, ynames, outnames, has_output,
                           return_by_column, return_initial,
                           return_time, return_output_with_y,
                           time_name) {
  bind_output <- has_output && return_output_with_y

  named <- FALSE
  if (has_output && !is.null(outnames)) {
    named <- return_output_with_y
    rownames(attr(ret, "output")) <- outnames
  }
  if (!is.null(ynames)) {
    named <- TRUE
    rownames(ret) <- ynames
  }

  if (return_time || bind_output) {
    at <- attributes(ret)
    if (return_time) {
      time <- matrix(if (return_initial) times else times[-1L], 1L,
                     dimnames = if (named) list(time_name, NULL) else NULL)
    } else {
      time <- NULL
    }

    ret <- rbind(if (return_time) time,
                 ret,
                 if (bind_output) at[["output"]],
                 deparse.level = 0L)
    if (bind_output) {
      at[["output"]] <- NULL
      has_output <- FALSE
    }
    ## This is a real pain, but we need to include any attributes set
    ## on the output by Cdopri; this is going to be "statistics" and
    ## "history", but it's always possible that additional attributes
    ## will be added later.
    for (x in setdiff(names(at), c("dim", "dimnames"))) {
      attr(ret, x) <- at[[x]]
    }
  }

  if (return_by_column) {
    ret <- t.default(ret)
    if (has_output) {
      attr(ret, "output") <- t.default(attr(ret, "output"))
    }
  }

  ret
}


find_function_address <- function(fun, dllname = "") {
  if (is.character(fun)) {
    fun <- getNativeSymbolInfo(fun, dllname)$address
  } else if (inherits(fun, "NativeSymbolInfo")) {
    fun <- fun$address
  } else if (!(inherits(fun, "externalptr") || is.function(fun))) {
    stop(sprintf("Invalid input for '%s'", deparse(substitute(fun))))
  }
  fun
}
