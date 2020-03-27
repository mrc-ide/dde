##' Solve a replicate set of difference (or recurrence) equation by
##' iterating it a number of times.  This is a wrapper around
##' \code{\link{difeq}} that does not (yet) do anything clever to
##' avoid many allocations.
##'
##' It is not currently possible to replicate over a set of parameters
##' at once yet; the same parameter set will be used for all
##' replications.
##'
##' The details of how replication is done here are all considered
##' implementation details and are up for change in the future - in
##' particular if the models are run in turn or simultaneously (and
##' the effect that has on the random number stream).  Logic around
##' naming output may change in future too; note that varying names in
##' the \code{y} here will have some unexpected behaviours.
##'
##' @title Solve difference equations repeatedly
##'
##' @param n Number of replicates.  It is an error to request zero
##'   replicates.
##'
##' @param y The initial state of the system.  Must be either a
##'   numeric vector or a \code{list} of numeric vectors.  If the
##'   latter, it must have length \code{n}.
##'
##' @inheritParams difeq
##'
##' @param as_array (Defunct) Logical, indicating if the output should
##'   be converted into an array.  If \code{TRUE} then \code{res[, ,
##'   i]} will contain the \code{i}'th replicate, if \code{FALSE} then
##'   \code{res[[i]]} does instead.  If both \code{as_array} and
##'   \code{restartable} are \code{TRUE}, then the attributes
##'   \code{ptr} and \code{restart_data} will be present as a
##'   \code{list} of restarting information for \code{difeq_continue},
##'   though using these is not yet supported.
##'
##' @export
##' @examples
##'
##' # Here is a really simple equation that does a random walk with
##' # steps that are normally distributed:
##' rhs <- function(i, y, p) y + runif(1)
##' y0 <- 1
##' t <- 0:10
##' p <- 5
##' dde::difeq_replicate(10, y0, t, rhs, p)
difeq_replicate <- function(n, y, steps, target, parms, ...,
                            n_out = 0L, n_history = 0L, grow_history = FALSE,
                            return_history = n_history > 0,
                            dllname = "", parms_are_real = TRUE,
                            ynames = NULL, outnames = NULL,
                            return_by_column = TRUE, return_initial = TRUE,
                            return_step = TRUE, return_output_with_y = TRUE,
                            restartable = FALSE, return_minimal = FALSE,
                            as_array = TRUE) {
  DOTS <- list(...)
  if (length(DOTS) > 0L) {
    stop("Invalid dot arguments!")
  }

  if (!as_array) {
    stop("as_array = FALSE no longer supported")
  }

  assert_positive_integer(n)

  if (is.list(y)) {
    assert_length(y, n)
    if (!all(lengths(y) == length(y[[1L]]))) {
      stop("All 'y' lengths must be the same")
    }
    nms <- lapply(y, names)
    if (length(unique(nms)) != 1) {
      stop("All 'y' names must be the same")
    }
    y <- matrix(unname(unlist(y)), length(y[[1L]]), length(y),
                dimnames = list(nms[[1L]], NULL))
  } else if (is.matrix(y)) {
    if (ncol(y) != n) {
      stop(sprintf("Expected '%d' columns in 'y'", n))
    }
  } else if (is.vector(y)) {
    y <- matrix(rep(y, n), length(y), n,
                dimnames = list(names(y), NULL))
  }
  if (storage.mode(y) != "double") {
    storage.mode(y) <- "double"
  }

  if (is.null(ynames) && !is.null(rownames(y))) {
    ynames <- rownames(y)
  } else {
    ynames <- check_ynames(y, ynames)
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
    parms <- list(target = target, parms = parms, rho = parent.frame())
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

  ret <- .Call(Cdifeq, y, as.integer(steps), target, parms,
               as.integer(n_out), parms_are_real,
               ## Return information:
               as.integer(n_history), grow_history, return_history,
               return_initial, restartable)

  has_output <- n_out > 0L
  prepare_output_replicate(ret, steps, ynames, outnames, has_output,
                           return_by_column, return_initial, return_step,
                           return_output_with_y,
                           "step")
}


prepare_output_replicate <- function(ret, times, ynames, outnames, has_output,
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
      attr(ret, "output") <- aperm(output, c(2, 1, 3))
    }
  }

  ret
}
