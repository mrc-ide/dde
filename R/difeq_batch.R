difeq_batch <- function(y, steps, target, parms, ...,
                        n_out = 0L, n_history = 0L, grow_history = FALSE,
                        return_history = n_history > 0,
                        dllname = "", parms_are_real = TRUE,
                        ynames = names(y), outnames = NULL,
                        return_by_column = TRUE, return_initial = TRUE,
                        return_step = TRUE, return_output_with_y = TRUE,
                        restartable = FALSE, return_minimal = FALSE) {
  if (!is.list(parms)) {
    ## TODO: could accept a matrix here and convert too
    stop("Expected a list of parameters")
  }
  if (is.list(y)) {
    len <- lengths(y)
    stopifnot(length(unique(len)) == 1L)
    y <- matrix(unlist(y), len, length(y))
  }

  ## TODO: for now start by replicating exactly the difeq starting
  ## point and we'l see where we get to.
  DOTS <- list(...)
  if (length(DOTS) > 0L) {
    stop("Invalid dot arguments!")
  }
  if (return_minimal) {
    return_by_column <- FALSE
    return_initial <- FALSE
    return_step <- FALSE
    return_output_with_y <- FALSE
  }

  target <- find_function_address(target, dllname)

  is_r_target <- is.function(target)

  if (is_r_target) {
    parms_are_real <- FALSE
    rho <- parent.frame()
    parms <- lapply(parms, function(p)
      list(target = target, parms = p, rho = rho))
    target <- NULL
    if (nzchar(dllname)) {
      stop("dllname must not be given when using an R function for 'target'")
    }
  }

  assert_size(n_history)
  assert_scalar_logical(grow_history)
  assert_scalar_logical(return_history)
  assert_scalar_logical(parms_are_real)
  assert_scalar_logical(return_by_column)
  assert_scalar_logical(return_initial)
  assert_scalar_logical(return_step)
  assert_scalar_logical(return_output_with_y)
  assert_scalar_logical(restartable)

  if (restartable) {
    stop("Restartable batch models not yet possible")
  }
  if (return_history) {
    stop("History-returning batch models not yet possible")
  }

  ## TODO: this will fail
  ynames <- check_ynames(y, ynames)

  assert_size(n_out)
  outnames <- check_outnames(n_out, outnames)

  if (steps[[1L]] < 0) {
    stop("steps must be positive")
  }
  if (length(steps) == 1L) {
    steps <- 0:steps
  }
  ## Avoid a class of issues
  if (n_history > 0 && n_history < 2L) {
    stop("If given, n_history must be at least 2")
  }

  ret <- .Call(Cdifeq_batch, y, as.integer(steps), target, parms,
               as.integer(n_out), parms_are_real,
               ## Return information:
               as.integer(n_history), grow_history, return_initial)
  ny <- if (is.matrix(y)) nrow(y) else length(y)
  dim(ret) <- c(ny, length(steps), length(parms))

  has_output <- n_out > 0L
  if (has_output) {
    dim(attr(ret, "output")) <- c(length(n_out), length(steps), length(parms))
  }

  ## Then we start with the repacking output.  This is actually pretty
  ## tedious to deal with
  prepare_output_batch(ret, steps, ynames, outnames, has_output,
                       return_by_column, return_initial, return_step,
                       return_output_with_y,
                       "step")
}


prepare_output_batch <- function(ret, times, ynames, outnames, has_output,
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
  }

  if (return_time || bind_output) {
    y <- ret
    at <- attributes(y)
    output <- attr(y, "output")
    n_out <- if (bind_output) nrow(output) else 0L

    ny <- dim(y)[[1L]]
    nt <- dim(y)[[2L]]
    np <- dim(y)[[3L]]

    ret <- array(NA_real_, c(as.integer(return_time) + ny + n_out, nt, np))
    ret[1, , ] <- rep(times, np)
    ret[seq_len(ny) + 1L, , ] <- y
    if (bind_output) {
      ret[seq_len(n_out) + ny + 1L, , ] <- output
      has_output <- FALSE
    }

    if (named) {
      rownames(ret) <- c(if (return_time) time_name,
                         ynames,
                         if (bind_output) outnames)
    }
  } else if (named) {
    rownames(ret) <- ynames
  }

  if (return_by_column) {
    ret <- aperm(ret, c(2, 1, 3))
    if (has_output) {
      attr(ret, "output") <- aperm(attr(ret, "output"), c(2, 1, 3))
    }
  }

  ret
}
