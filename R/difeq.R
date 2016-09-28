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
##' @param return_step Logical, indicating if a row (or column if
##'   \code{by_column} is \code{TRUE}) representing step is included.
##'
##' @param restartable Logical, indicating if the problem should be
##'   restartable.  If \code{TRUE}, then the return value of a
##'   simulation can be passed to \code{difeq_restart} to continue the
##'   simulation after arbitrary changes to the state or the
##'   parameters.  Note that this is really only useful for delay
##'   difference equations where you want to keep the history but make
##'   changes to the parameters or to the state vector while keeping
##'   the history of the problem so far.
##'
##' @inheritParams dopri
##' @export
difeq <- function(y, steps, target, parms, ...,
                  n_out = 0L, n_history = 0L, grow_history = FALSE,
                  return_history = n_history > 0,
                  dllname = "", parms_are_real = TRUE,
                  ynames = TRUE, outnames = NULL,
                  by_column = FALSE, return_initial = FALSE,
                  return_step = FALSE, return_output_with_y = FALSE,
                  restartable = FALSE, deSolve_compatible = FALSE) {
  DOTS <- list(...)
  if (length(DOTS) > 0L) {
    stop("Invalid dot arguments!")
  }
  if (deSolve_compatible) {
    by_column <- TRUE
    return_initial <- TRUE
    return_step <- TRUE
    return_output_with_y <- TRUE
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
  assert_scalar_logical(grow_history)
  assert_scalar_logical(return_history)
  assert_scalar_logical(parms_are_real)
  assert_scalar_logical(by_column)
  assert_scalar_logical(return_initial)
  assert_scalar_logical(return_step)
  assert_scalar_logical(return_output_with_y)
  assert_scalar_logical(restartable)

  ynames <- check_ynames(y, ynames, deSolve_compatible)

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
  ret <- prepare_output(ret, steps, ynames, outnames, has_output,
                        by_column, return_initial, return_step,
                        return_output_with_y,
                        "step")
  if (restartable) {
    ret <- prepare_difeq_restart(ret, parms, parms_are_real, ynames, outnames,
                                 has_output, return_history, by_column,
                                 return_initial, return_step,
                                 return_output_with_y)
  }
  ret
}

##' @export
##' @rdname difeq
##' @param obj An object to continue from; this must be the results of
##'   running a simulation with the option \code{restartable =
##'   TRUE}.  Note that continuing a problem moves the pointer along
##'   in time (unless \code{copy = TRUE}, and that the incoming time
##'   (\code{times[[1]]}) must equal the previous time \emph{exactly}.
##'
##' @param copy Logical, indicating if the pointer should be copied
##'   before continuing.  If \code{TRUE}, this is non-destructive with
##'   respect to the data in the original pointer so the problem can
##'   be restarted multiple times.  By default this is \code{FALSE}
##'   becaue there is a (potentially very small) cost to this
##'   operation.
difeq_continue <- function(obj, steps, y = NULL, ...,
                           copy = FALSE,
                           parms = NULL,
                           return_history = NULL,
                           by_column = NULL, return_initial = NULL,
                           return_time = NULL, return_output_with_y = NULL,
                           restartable = NULL) {
  DOTS <- list(...)
  if (length(DOTS) > 0L) {
    stop("Invalid dot arguments!")
  }

  ptr <- attr(obj, "ptr", exact = TRUE)
  dat <- attr(obj, "restart_data", exact = TRUE)
  if (copy) {
    ptr <- .Call(Cdopri_copy, ptr)
  }

  if (is.null(parms)) {
    parms <- dat$parms
  }

  ## Process any given options, falling back on the previous values
  return_history <- logopt(return_history, dat$return_history)
  by_column <- logopt(by_column, dat$by_column)
  return_initial <- logopt(return_initial, dat$return_initial)
  return_time <- logopt(return_time, dat$return_time)
  return_output_with_y <- logopt(return_output_with_y, dat$return_output_with_y)
  restartable <- logopt(restartable, TRUE)

  ret <- .Call(Cdifeq_continue, ptr, y, as.integer(steps),
               parms, dat$parms_are_real,
               return_history, return_initial, return_statistics, restartable)

  ret <- prepare_output(ret, step, dat$ynames, dat$outnames, dat$has_output,
                        by_column, return_initial, return_time,
                        return_output_with_y, "step")
  if (restartable) {
    ret <- prepare_difeq_restart(ret, parms, parms_are_real, ynames, outnames,
                                 has_output, return_history, by_column,
                                 return_initial, return_step,
                                 return_output_with_y)
  }
  ret
}

##' @export
##' @rdname difeq
##'
##' @param step The step to access (not that this is not an offset,
##'   but the actual step; within your target function you'd write
##'   things like \code{yprev(step - 1)} to get the previous step.
##'
##' @param i index within the state vector \code{y} to return.  The
##'   index here is R-style base-1 indexing, so pass \code{1} in to
##'   access the first element.  This can be left \code{NULL} to
##'   return all the elements or a vector longer than one.
yprev <- function(step, i = NULL) {
  .Call(Cyprev, step, i)
}

prepare_difeq_restart <- function(ret, parms, parms_are_real, ynames, outnames,
                                  has_output, return_history, by_column,
                                  return_initial, return_step,
                                  return_output_with_y) {
  restart_data <- list(parms = parms, parms_are_real = parms_are_real,
                       ynames = ynames, outnames = outnames,
                       has_output = has_output,
                       return_history = return_history,
                       by_column = by_column,
                       return_initial = return_initial,
                       return_step = return_step,
                       return_output_with_y = return_output_with_y)
  attr(ret, "restart_data") <- restart_data
  ret
}
