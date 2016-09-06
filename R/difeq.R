##' Solve a difference (or recurrence) equation by iterating it a
##' number of times.
##'
##' @title Solve difference equation
##'
##' @param y The initial state of the system.  Must be a numeric vector
##'
##' @param steps Either a vector of steps to solve the system at or a
##'   positive integer.  If a vector of steps, then the \emph{first}
##'   step is taken as step zero, and the solution will be recorded at
##'   every other step in the vector.  So to step a system from time
##'   zero to times 1, 2, 3, ... use 0:n.
##'
##' @param target The target function to advance.  This can either be
##'   an R function taking arguments \code{n, i, t, y, parms} or be a
##'   scalar character with the name of a compiled function with
##'   arguments \code{size_t n, size_t step, double time, const double
##'   *y, double *dydt, size_t n_out, double *output void *data}.
##
##' @param parms Parameters to pass through to the difference function
##'
##' @param n_out Number of "output" variables (not differential
##'   equation variables) to compute via the routine \code{output}.
##'
##' @param t0 The initial time, if your problem has a continuous time
##'   basis too (e.g., if each step corresponds to 1/365 of a year)
##'
##' @param dt The continuous time step size
##'
##' @param n_out The number of output variables (in addition to the
##'   difference equation variables).  If given, then an R function
##'   must return an \emph{attribute} \code{output} with the output
##'   variables.
##'
##' @param n_history The number of iterations of history to save
##'   during the simulation.  By default, no history is saved.
##'
##' @param return_history Logical indicating if history is to be
##'   returned.  By default, history is returned if \code{n_history}
##'   is nonzero.
##'
##' @param dllname Name of the shared library (without extension) to
##'   find the function \code{func} in the case where \code{func}
##'   refers to compiled function.
##'
##' @inheritParams dopri
##' @export
difeq <- function(y, steps, target, parms, ...,
                  t0 = 0, dt = 1,
                  n_out = 0L, n_history = 0L, return_history = n_history > 0,
                  dllname = "", parms_are_real = TRUE,
                  ynames = TRUE, outnames = NULL,
                  by_column = FALSE, return_initial = FALSE,
                  return_step = FALSE, deSolve_compatible = FALSE) {
  DOTS <- list(...)
  if (length(DOTS) > 0L) {
    stop("Invalid dot arguments!")
  }
  if (deSolve_compatible) {
    by_column <- TRUE
    return_initial <- TRUE
    return_step <- TRUE
  }

  target <- find_function_address(target, dllname)

  is_r_target <- is.function(target)

  if (is_r_target) {
    parms_are_real <- FALSE
    parms <- list(target = target, parms = parms, rho = parent.frame())
    target <- find_function_address("difeq_r_harness", "dde")
    if (nzchar(dllname)) {
      stop("dllname must not be given when using an R function for 'target'")
    }
  }

  assert_size(n_history)
  assert_scalar_logical(return_history)
  assert_scalar_logical(parms_are_real)
  assert_scalar_logical(by_column)
  assert_scalar_logical(return_initial)
  assert_scalar_numeric(t0)
  assert_scalar_numeric(dt)

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

  if (return_step && !is.null(ynames)) {
    ynames <- c("step", ynames)
  }

  assert_size(n_out)
  if (!is.null(outnames)) {
    if (is.character(outnames)) {
      if (length(outnames) != n_out) {
        stop("outnames must have length n_out")
      }
    } else {
      stop("Invalid value for outnames")
    }
  }

  if (n_out > 0L && is_r_target) {
    stop("Not supported yet")
  }

  if (any(steps < 0)) {
    stop("steps must be positive")
  }
  if (length(steps) == 1L) {
    steps <- 0:steps
  }

  ret <- .Call(Cdifeq, y, as.integer(steps), target, parms,
               as.integer(n_out), as.numeric(t0), as.numeric(dt),
               parms_are_real,
               ## Return information:
               as.integer(n_history), return_history,
               return_initial)

  if (return_step) {
    at <- attributes(ret)
    ret <- rbind(if (return_initial) steps else steps[-1],
                 ret, deparse.level = 0)
    ## This is a real pain, but we need to include any attributes set
    ## on the output by Cdopri; this is going to be "statistics" and
    ## "history", but it's always possible that additional attributes
    ## will be added.
    for (x in setdiff(names(at), "dim")) {
      attr(ret, x) <- at[[x]]
    }
  }

  if (!is.null(ynames)) {
    rownames(ret) <- ynames
  }
  if (n_out > 0L && !is.null(outnames)) {
    rownames(attr(ret, "output")) <- outnames
  }

  if (by_column) {
    ret <- t.default(ret)
    if (n_out > 0L) {
      attr(ret, "output") <- t.default(attr(ret, "output"))
    }
  }

  ret
}
