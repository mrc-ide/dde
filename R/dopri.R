##' Integrate an ODE or DDE with dopri.
##'
##' Like \code{deSolve::lsoda}, this function has \emph{many}
##' arguments.  This is far from ideal, and I would welcome any
##' approach for simplifying it a bit.
##'
##' The options \code{return_by_column}, \code{return_initial},
##' \code{return_time}, \code{return_output_with_y} exist because
##' these options all carry out modifications of the data at the end
##' of solving the ODE and this can incur a small but measurable cost.
##' When solving an ODE repeatedly (e.g., in the context of an MCMC or
##' optimisation) it may be useful to do as little as possible.  For
##' simple problems this can save around 5-10\% of the total
##' computational time (especially the transpose).  The shorthand
##' option \code{return_minimal} will set all to \code{FALSE} when
##' used.
##'
##' @title Integrate ODE/DDE with dopri
##'
##' @param y Initial conditions for the integration
##'
##' @param times Times where output is needed.  Unlike \code{deSolve}
##'   we won't actually stop at these times, but instead interpolate
##'   back to get the result.
##'
##' @param func Function to integrate.  Can be an R function of
##'   arguments \code{t, y, parms}, returning a numeric vector, or it
##'   can be the name or address of a C function with arguments
##'   \code{size_t n, double t, const double *y, double *dydt, void *data}.
##'
##' @param parms Parameters to pass through to the derivatives.
##'
##' @param ... Dummy arguments - nothing is allowed here, but this
##'   means that all further arguments \emph{must} be specified by
##'   name (not order) so I can easily reorder them later on.
##'
##' @param n_out Number of "output" variables (not differential
##'   equation variables) to compute via the routine \code{output}.
##'
##' @param output The output routine; either an R function taking
##'   arguments \code{t, y, parms} or the name/address of a C function
##'   taking arguments \code{size_t n, double t, const double *y,
##'   size_t n_out, double *out, void *data}.
##'
##' @param rtol The per-step relative tolerance.  The total accuracy
##'   will be less than this.
##'
##' @param atol The per-step absolute tolerance.
##'
##' @param step_size_min The minimum step size.  The actual minimum
##'   used will be the largest of the absolute value of this
##'   \code{step_size_min} or \code{.Machine$double.eps}.  If the
##'   integration attempts to make a step smaller than this, it will
##'   throw an error by default, stopping the integration (note that
##'   this differs from the treatment of \code{hmin} in
##'   \code{deSolve::lsoda}). See \code{allow_step_size_min} to change
##'   this behaviour.
##'
##' @param step_size_max The largest step size.  By default there is
##'   no maximum step size (Inf) so the solver can take as large a
##'   step as it wants to.  If you have short-lived fluctuations in
##'   your rhs that the solver may skip over by accident, then specify
##'   a smaller maximum step size here (or use \code{tcrit} below).
##'
##' @param step_size_initial The initial step size.  By default the
##'   integrator will guess the step size automatically, but one can
##'   be given here instead.
##'
##' @param step_max_n The maximum number of steps allowed.  If the
##'   solver takes more steps than this it will throw an error.  Note
##'   the number of evaluations of \code{func} will be about 6 times
##'   the number of steps (or 11 times if using \code{method =
##'   "dopri853"}).
##'
##' @param step_size_min_allow Logical, indicating if when a step size
##'   is driven down to \code{step_size_min} we should allow it to
##'   proceed. This is the behaviour in of \code{hmin} in
##'   \code{deSolve::lsoda}.
##'
##' @param tcrit An optional vector of critical times that the solver
##'   must stop at (rather than interpolating over).  This can include
##'   an end time that we can't go past, or points within the
##'   integration that must be stopped at exactly (for example cases
##'   where the derivatives change abruptly).  Note that this differs
##'   from the interpretation of this parameter in deSolve; there
##'   \code{tcrit} is a single time that integration may not go past
##'   -- with dde we never go past the final time, and this is just
##'   for times that fall \emph{within} the range of times in
##'   \code{times}.
##'
##' @param event_time Vector of times to fire events listed in
##'   \code{event_function} at
##'
##' @param event_function Function to fire at events.  For R models
##'   (\code{func} is an R function and \code{dllname} is empty), this
##'   must be either a single R function (same function for all
##'   events) or a \code{list} of R functions.  For C models, this
##'   must be a singe C function (same requirements as \code{func} or
##'   \code{output} or a list/vector of these as appropriate).
##'
##' @param method The integration method to use, as a string.  The
##'   supported methods are \code{"dopri5"} (5th order method with 4th
##'   order dense output) and \code{"dopri853"} (8th order method with
##'   7th order output and embedded 5th and 3rd order schemes).
##'   Alternatively, use the functions \code{dopri5} or
##'   \code{dopri853} which simply sets this argument.
##'
##' @param stiff_check How often to check that the problem has become
##'   stiff.  If zero, then the problem is never checked, and if
##'   positive then the problem is checked every \code{stiff_check}
##'   accepted steps.  The actual check is based off the algorithm in
##'   Hairer's implementation of the solvers and may be overly strict,
##'   especially for delay equations with the 853 method (in my
##'   limited experience with it).
##'
##' @param verbose Be verbose, and print information about each step.
##'   This may be useful for learning about models that misbehave.
##'   Valid values are \code{TRUE} (enable debugging) or \code{FALSE}
##'   (disable debugging) or use one of \code{dopri:::VERBOSE_QUIET},
##'   \code{dopri:::VERBOSE_STEP} or \code{VERBOSE:::VERBOSE_EVAL}.
##'   If an R function is provided as the argument \code{callback}
##'   then this function will also be called at each step or
##'   evaluation (see below for details).
##'
##' @param callback Callback function that can be used to make verbose
##'   output more useful.  This can be used to return more information
##'   about the evaluation as it proceeds, generally as information
##'   printed to the screen.  The function must accept arguments
##'   \code{t}, \code{y} and \code{dydt}.  See Details for further
##'   information.
##'
##' @param n_history Number of history points to retain.  This needs
##'   to be greater than zero for delay differential equations to
##'   work.  Alternatively, this may be greater than zero to return
##'   model outputs that can be inspected later.
##'
##' @param grow_history Logical indicating if history should be grown
##'   during the simulation.  If \code{FALSE} (the default) then when
##'   history is used it is overwritten as needed (so only the most
##'   recent \code{n_history} elements are saved.  This may require
##'   some tuning so that you have enough history to run your
##'   simulation (i.e. to the longest delay) or an error will be
##'   thrown when it underflows.  The required history length will
##'   vary with your delay sizes and with the timestep for dopri.  If
##'   \code{TRUE}, then history will grow as the buffer is exhausted.
##'   The growth is geometric, so every time it reaches the end of the
##'   buffer it will increase by a factor of about 1.6 (see the
##'   \code{ring} documentation).  This may consume more memory than
##'   necessary, but may be useful where you don't want to care about
##'   picking the history length carefully.
##'
##' @param return_history Logical indicating if history should be
##'   returned alongside the output or discarded.  By default, history
##'   is retained if \code{n_history} is greater than 0, but that
##'   might change (and may not be desirable unless you plan on
##'   actually using it).
##'
##' @param dllname Name of the shared library (without extension) to
##'   find the function \code{func} (and \code{output} if given) in
##'   the case where \code{func} refers to compiled function.
##'
##' @param parms_are_real Logical, indicating if \code{parms} should
##'   be treated as vector of doubles by \code{func} (when it is a
##'   compiled function).  If \code{TRUE} (the default), then
##'   \code{REAL(parms)}, which is \code{double*} is passed through.
##'   If \code{FALSE} then if \code{params} is an externalptr type
##'   (\code{EXTPTRSXP}) we pass through the result of
##'   \code{R_ExternalPtrAddr}, otherwise we pass \code{params}
##'   through unmodified as a \code{SEXP}.  In the last case, in your
##'   target function you will need to include \code{<Rinternals.h>},
##'   cast to \code{SEXP} and then pull it apart using the R API (or
##'   Rcpp).
##'
##' @param ynames Logical, indicating if the output should be named
##'   following the names of the input vector \code{y}.
##'   Alternatively, if \code{ynames} is a character vector of the
##'   same length as \code{y}, these will be used as the output names.
##'
##' @param outnames An optional character vector, used when
##'   \code{n_out} is greater than 0, to name the model output matrix.
##'
##' @param return_by_column Logical, indicating if the output should be
##'   returned organised by column (rather than row).  This incurs a
##'   slight cost for transposing the matrices.  If you can work with
##'   matrices that are transposed relative to \code{deSolve}, then
##'   set this to \code{FALSE}.
##'
##' @param return_initial Logical, indicating if the output should
##'   include the initial conditions.  Specifying \code{FALSE} avoids
##'   binding this onto the output.
##'
##' @param return_time Logical, indicating if a row (or column if
##'   \code{return_by_column} is \code{TRUE}) representing time is included.
##'   If \code{FALSE}, this is not added.
##'
##' @param return_output_with_y Logical, indicating if the output
##'   should be bound together with the returned matrix \code{y} (as
##'   it is with \code{deSolve}).  If \code{FALSE}, then output will
##'   be returned as the attribute \code{output}.
##'
##' @param return_statistics Logical, indicating if statistics about
##'   the run should be included.  If \code{TRUE}, then an integer
##'   vector containing the number of target evaluations, steps,
##'   accepted steps and rejected steps is returned (the vector is
##'   named).
##'
##' @param return_minimal Shorthand option - if set to \code{TRUE}
##'   then it sets all of \code{return_by_column},
##'   \code{return_initial}, \code{return_time},
##'   \code{return_output_with_y} to \code{FALSE}
##'
##' @param restartable Logical, indicating if the problem should be
##'   restartable.  If \code{TRUE}, then the return value of an
##'   integration can be passed to \code{dopri_restart} to continue
##'   the integration after arbitrary changes to the state or the
##'   parameters.  Note that when using delay models, the integrator
##'   is fairly naive about how abrupt changes in the state space are
##'   dealt with, and may perform very badly with \code{method =
##'   "dopri853"} which assumes a fairly smooth problem.  Note that
##'   this is really only useful for delay differential equations
##'   where you want to keep the history but make changes to the
##'   parameters or to the state vector while keeping the history of
##'   the problem so far.
##'
##' @return At present the return value is transposed relative to
##'   deSolve.  This might change in future.
##'
##' @export
##'
##' @seealso \code{\link{dopri_interpolate}} which can be used to
##'   efficiently sample from output of \code{dopri}, and the package
##'   vignette which shows in more detail how to solve delay
##'   differential equations and to use compiled objective functions.
##'
##' @section Verbose output and callbacks:
##'
##' Debugging a failed integration can be difficult, but \code{dopri}
##'   provides a couple of tools to get more information about where a
##'   failure might have occurred.  Most simply, one can pass
##'   \code{verbose = TRUE} which will print information about the
##'   time and the step size at each point just before the step is
##'   stated.  Passing in \code{verbose = dde:::VERBOSE_EVAL} will
##'   print information just before every evaluation of the target
##'   function (there are several evaluations per step).
##'
##' However, this does not provide information about the state just
##'   before failure.  To get that, one must provide a \code{callback}
##'   function - this is an R function that will be called just before
##'   a step or evaluation (based on the value of the \code{verbose}
##'   argument) in place of the default print.  Define a callback
##'   function with arguments \code{t}, \code{h} and \code{y} where
##'   \code{t} is the time (beginning of a step or location of an
##'   evaluation), \code{h} is the step size (or \code{NA} for an
##'   evaluation) and \code{y} is the state at the point of the step
##'   or evaluation.  Your callback function can do anything - you can
##'   print to the screen (using \code{cat} or \code{message}), you
##'   can store results using a closure and \code{<<-} or you could
##'   conditionally use a \code{browser()} call to debug
##'   interactively.  However, it is not possible for the callback to
##'   affect the solution short of throwing an error and interrupting
##'   it.  See the Examples for an example of use.
##'
##' @examples
##'
##' # The lorenz attractor:
##' lorenz <- function(t, y, p) {
##'   sigma <- p[[1L]]
##'   R <- p[[2L]]
##'   b <- p[[3L]]
##'   c(sigma * (y[[2L]] - y[[1L]]),
##'     R * y[[1L]] - y[[2L]] - y[[1L]] * y[[3L]],
##'     -b * y[[3L]] + y[[1L]] * y[[2L]])
##' }
##'
##' p <- c(10, 28, 8 / 3)
##' y0 <- c(10, 1, 1)
##'
##' tt <- seq(0, 100, length.out = 40000)
##' y <- dde::dopri(y0, tt, lorenz, p, return_time = FALSE)
##' plot(y[, c(1, 3)], type = "l", lwd = 0.5, col = "#00000066")
##'
##' # If we want to print progress as the integration progresses we can
##' # use the verbose argument:
##' y <- dde::dopri(y0, c(0, 0.1), lorenz, p, verbose = TRUE)
##'
##' # Or print the y values too using a callback:
##' callback <- function(t, h, y) {
##'   message(sprintf("t: %f, h: %e, y: [%s]", t, h,
##'                   paste(format(y, 5), collapse = ", ")))
##' }
##' y <- dde::dopri(y0, c(0, 0.1), lorenz, p, verbose = TRUE,
##'                 callback = callback)
dopri <- function(y, times, func, parms, ...,
                  n_out = 0L, output = NULL,
                  rtol = 1e-6, atol = 1e-6,
                  step_size_min = 0, step_size_max = Inf,
                  step_size_initial = 0, step_max_n = 100000L,
                  step_size_min_allow = FALSE,
                  tcrit = NULL, event_time = NULL, event_function = NULL,
                  method = "dopri5",
                  stiff_check = 0,
                  verbose = FALSE, callback = NULL,
                  n_history = 0, grow_history = FALSE,
                  return_history = n_history > 0, dllname = "",
                  parms_are_real = TRUE,
                  ynames = names(y), outnames = NULL,
                  return_by_column = TRUE, return_initial = TRUE,
                  return_time = TRUE, return_output_with_y = TRUE,
                  return_statistics = FALSE, restartable = FALSE,
                  return_minimal = FALSE) {
  ## TODO: include "deSolve" mode where we do the transpose, add the
  ## time column too?
  DOTS <- list(...)
  if (length(DOTS) > 0L) {
    stop("Invalid dot arguments!")
  }
  if (return_minimal) {
    return_by_column <- FALSE
    return_initial <- FALSE
    return_time <- FALSE
    return_output_with_y <- FALSE
    ynames <- NULL
  }

  func <- find_function_address(func, dllname)

  is_r_target <- is.function(func)

  if (is_r_target) {
    parms_are_real <- FALSE
    parms <- list(func = func, parms = parms, rho = parent.frame())
    func <- NULL
    if (nzchar(dllname)) {
      stop("dllname must not be given when using an R function for 'func'")
    }
  }

  use_853 <- match_value(method, dopri_methods()) == "dopri853"

  assert_scalar(rtol)
  assert_scalar(atol)
  assert_scalar(step_size_min)
  assert_scalar(step_size_max)
  assert_scalar(step_size_initial)
  assert_size(step_max_n)
  assert_scalar_logical(step_size_min_allow)
  assert_size(n_history)

  assert_scalar_logical(grow_history)
  assert_scalar_logical(return_history)
  assert_scalar_logical(parms_are_real)
  assert_scalar_logical(return_by_column)
  assert_scalar_logical(return_initial)
  assert_scalar_logical(return_statistics)
  assert_scalar_logical(return_time)
  assert_scalar_logical(return_output_with_y)
  assert_scalar_logical(restartable)
  assert_size(stiff_check)

  verbose <- dopri_verbose(verbose)
  callback <- dopri_callback(callback)

  ynames <- check_ynames(y, ynames)

  assert_size(n_out)
  outnames <- check_outnames(n_out, outnames)

  if (n_out > 0L) {
    output <- find_function_address(output, dllname)
    ## NOTE: The same-typedness of output/func is not really
    ## necessary, but I think it's simplest to think about if we
    ## enforce it.
    if (is_r_target) {
      if (!is.function(output)) {
        stop("output must be an R function")
      }
      parms[[DOPRI_IDX_OUTPUT]] <- output
      output <- NULL
    } else {
      if (is.function(output)) {
        stop("output must be a compiled function (name or address)")
      }
    }
  } else if (!is.null(output)) {
    stop("If 'output' is given, then n_out must be specified")
  }

  dat <- check_events(event_time, event_function, tcrit, dllname)
  tcrit <- dat$tcrit
  ## This is needed to allow things like `tcrit = 1:5` (which is int)
  if (!is.null(tcrit)) {
    assert_numeric(tcrit)
    tcrit <- as.numeric(tcrit)
  }

  is_event <- dat$is_event
  event <- dat$event

  if (!is.null(event) && is_r_target) {
    parms[[DOPRI_IDX_EVENT]] <- dat$event_r
  }

  ret <- .Call(Cdopri, y, as.numeric(times), func, parms,
               as.integer(n_out), output, parms_are_real,
               ## Tolerance:
               rtol, atol,
               ## Step control:
               step_size_min, step_size_max,
               step_size_initial, as.integer(step_max_n),
               step_size_min_allow,
               ## Critical and events
               as.numeric(tcrit), is_event, event,
               ## Other:
               use_853, as.integer(stiff_check), verbose, callback,
               ## Return information:
               as.integer(n_history), grow_history, return_history,
               return_initial, return_statistics, restartable)

  has_output <- n_out > 0L
  ret <- prepare_output(ret, times, ynames, outnames, has_output,
                        return_by_column, return_initial, return_time,
                        return_output_with_y,
                        "time")
  if (restartable) {
    restart_data <- list(parms = parms, parms_are_real = parms_are_real,
                         ynames = ynames, outnames = outnames,
                         has_output = has_output,
                         tcrit = tcrit,
                         return_history = return_history,
                         return_by_column = return_by_column,
                         return_initial = return_initial,
                         return_statistics = return_statistics,
                         return_time = return_time,
                         return_output_with_y = return_output_with_y)
    attr(ret, "restart_data") <- restart_data
  }
  ret
}


##' @export
##' @rdname dopri
dopri5 <- function(y, times, func, parms, ...) {
  dopri(y, times, func, parms, ..., method = "dopri5")
}


##' @export
##' @rdname dopri
dopri853 <- function(y, times, func, parms, ...) {
  dopri(y, times, func, parms, ..., method = "dopri853")
}


##' @export
##' @rdname dopri
##' @param obj An object to continue from; this must be the results of
##'   running an integration with the option \code{restartable =
##'   TRUE}.  Note that continuing a problem moves the pointer along
##'   in time (unless \code{copy = TRUE}, and that the incoming time
##'   (\code{times[[1]]}) must equal the previous time \emph{exactly}.
##'
##' @param copy Logical, indicating if the pointer should be copied
##'   before continuing.  If \code{TRUE}, this is non-destructive with
##'   respect to the data in the original pointer so the problem can
##'   be restarted multiple times.  By default this is \code{FALSE}
##'   because there is a (potentially very small) cost to this
##'   operation.
dopri_continue <- function(obj, times, y = NULL, ...,
                           copy = FALSE,
                           parms = NULL,
                           tcrit = NULL, return_history = NULL,
                           return_by_column = NULL, return_initial = NULL,
                           return_statistics = NULL, return_time = NULL,
                           return_output_with_y = NULL, restartable = NULL) {
  DOTS <- list(...)
  if (length(DOTS) > 0L) {
    stop("Invalid dot arguments!")
  }

  ptr <- attr(obj, "ptr", exact = TRUE)
  dat <- attr(obj, "restart_data", exact = TRUE)
  if (copy) {
    ptr <- .Call(Cdopri_copy, ptr)
  }

  if (is.null(tcrit)) {
    tcrit <- dat$tcrit
  }
  if (is.null(parms)) {
    parms <- dat$parms
  }

  ## Process any given options, falling back on the previous values
  return_history <- logopt(return_history, dat$return_history)
  return_by_column <- logopt(return_by_column, dat$return_by_column)
  return_initial <- logopt(return_initial, dat$return_initial)
  return_statistics <- logopt(return_statistics, dat$return_statistics)
  return_time <- logopt(return_time, dat$return_time)
  return_output_with_y <- logopt(return_output_with_y, dat$return_output_with_y)
  restartable <- logopt(restartable, TRUE)

  ret <- .Call(Cdopri_continue, ptr, y, as.numeric(times),
               parms, dat$parms_are_real, tcrit,
               return_history, return_initial, return_statistics, restartable)

  ## TODO: make this work as restartable
  prepare_output(ret, times, dat$ynames, dat$outnames, dat$has_output,
                 return_by_column, return_initial, return_time,
                 return_output_with_y, "time")
}


##' Interpolate the Dormand-Prince output after an integration.  This
##' only interpolates the core integration variables and not any
##' additional output variables.
##'
##' This decouples the integration of the equations and the generation
##' of output; it is not necessary for use of the package, but may
##' come in useful where you need to do (for example) root finding on
##' the time course of a problem, or generate minimal output in some
##' cases and interrogate the solution more deeply in others.  See the
##' examples and the package vignette for a full worked example.
##'
##' @title Interpolate Dormand-Prince output
##' @param h The interpolation history.  This can be the output
##'   running \code{dopri} with \code{return_history = TRUE}, or the
##'   history attribute of this object (retrievable with
##'   \code{attr(res, "history")}).
##' @param t The times at which interpolated output is required.
##'   These times must fall within the included history (i.e., the
##'   times that the original simulation was run) or an error will be
##'   thrown.
##' @export
##' @author Rich FitzJohn
##'
##' @examples
##' # Here is the Lorenz attractor implemented as an R function
##' lorenz <- function(t, y, p) {
##'   sigma <- p[[1L]]
##'   R <- p[[2L]]
##'   b <- p[[3L]]
##'   c(sigma * (y[[2L]] - y[[1L]]),
##'     R * y[[1L]] - y[[2L]] - y[[1L]] * y[[3L]],
##'     -b * y[[3L]] + y[[1L]] * y[[2L]])
##' }
##'
##' # Standard parameters and a reasonable starting point:
##' p <- c(10, 28, 8 / 3)
##' y0 <- c(10, 1, 1)
##'
##' # Run the integration for times [0, 50] and return minimal output,
##' # but *do* record and return history.
##' y <- dopri(y0, c(0, 50), lorenz, p,
##'            n_history = 5000, return_history = TRUE,
##'            return_time = FALSE, return_initial = FALSE,
##'            return_by_column = FALSE)
##'
##' # Very little output is returned (just 3 numbers being the final
##' # state of the system), but the "history" attribute is fairly
##' # large matrix of history information.  It is not printed though
##' # as its contents should not be relied on.  What does matter is
##' # the range of supported times printed (i.e., [0, 50]) and the
##' # number of entries (~2000).
##' y
##'
##' # Generate an interpolated set of variables using this; first for
##' # 1000 steps over the full range:
##' tt <- seq(0, 50, length.out = 1000)
##' yy <- dopri_interpolate(y, tt)
##' plot(yy[, c(1, 3)], type = "l")
##'
##' # Then for 50000
##' tt <- seq(0, 50, length.out = 50000)
##' yy <- dopri_interpolate(y, tt)
##' plot(yy[, c(1, 3)], type = "l")
dopri_interpolate <- function(h, t) {
  if (!inherits(h, "dopri_history")) {
    h <- attr(h, "history", exact = TRUE)
  }

  nd <- attr(h, "n") # number of equations
  if (is.null(nd)) {
    stop("Corrupt history object: 'n' is missing")
  }
  nh <- ncol(h)      # number of history entries
  it <- nrow(h) - 1L # time index
  ih <- nrow(h)      # step size index
  order <- (nrow(h) - 2) / nd # order of integration
  if (!(order == 5 || order == 8)) {
    ## This one should really never be triggered but acts as a
    ## safeguard against real weirdness with subsetting.
    stop("Corrupt history object: incorrect number of rows")
  }

  tr <- c(h[it, 1L], h[it, nh] + h[ih, nh])
  if (min(t) < tr[[1L]] || max(t) > tr[[2L]]) {
    stop(sprintf("Time falls outside of range of known history [%s, %s]",
                 tr[[1L]], tr[[2L]]))
  }

  idx <- findInterval(t, h[it, ])
  theta  <- (t - h[it, idx]) / h[ih, idx]
  theta1 <- 1.0 - theta

  ret <- matrix(NA_real_, length(t), nd)

  ## Next, we need to put history into the right shape.
  history <- array(h[seq_len(order * nd), ], c(nd, order, nh))
  ## Then pulling apart into the 5 coefficients
  h1 <- history[, 1L, ]
  h2 <- history[, 2L, ]
  h3 <- history[, 3L, ]
  h4 <- history[, 4L, ]
  h5 <- history[, 5L, ]
  if (order == 5) {
    for (i in seq_len(nd)) {
      ret[, i] <- h1[i, idx] + theta *
        (h2[i, idx] + theta1 *
         (h3[i, idx] + theta *
          (h4[i, idx] + theta1 *
           h5[i, idx])))
    }
  } else { # order == 8
    h6 <- history[, 6L, ]
    h7 <- history[, 7L, ]
    h8 <- history[, 8L, ]
    for (i in seq_len(nd)) {
      tmp <- h5[i, idx] + theta *
        (h6[i, idx] + theta1 *
         (h7[i, idx] + theta *
          h8[i, idx]))
      ret[, i] <- h1[i, idx] + theta *
        (h2[i, idx] + theta1 *
         (h3[i, idx] + theta *
          (h4[i, idx] + theta1 *
           tmp)))
    }
  }

  ret
}


##' @export
print.dopri_history <- function(x, ...) {
  nd <- attr(x, "n") # number of equations
  nh <- ncol(x)
  it <- nrow(x) - 1L
  ih <- nrow(x)
  order <- (nrow(x) - 2) / nd
  cat("<dopri_history>\n")
  cat(sprintf("  - equations: %d\n", nd))
  cat(sprintf("  - time: [%s, %s]\n", x[it, 1L], x[it, nh] + x[ih, nh]))
  cat(sprintf("  - entries: %d\n", nh))
  cat(sprintf("  - order: %d\n", order))
}


##' @export
##' @rdname dopri
##'
##' @param t The time to access (not that this is not an offset,
##'   but the actual time; within your target function you'd write
##'   things like \code{tlag(t - 1)} to get 1 time unit ago.
##'
##' @param i index within the state vector \code{y} to return.  The
##'   index here is R-style base-1 indexing, so pass \code{1} in to
##'   access the first element.  This can be left \code{NULL} to
##'   return all the elements or a vector longer than one.
ylag <- function(t, i = NULL) {
  .Call(Cylag, t, i)
}


dopri_methods <- function() {
  c("dopri5", "dopri853")
}


dopri_verbose <- function(verbose) {
  assert_scalar(verbose)
  if (is.logical(verbose)) {
    verbose <- if (verbose) VERBOSE_STEP else VERBOSE_QUIET
  } else {
    assert_integer(verbose)
    valid <- c(VERBOSE_QUIET, VERBOSE_STEP, VERBOSE_EVAL)
    if (!any(verbose == valid)) {
      stop("Invalid value for verbose")
    }
    verbose <- as.integer(verbose)
  }
  verbose
}


dopri_callback <- function(callback) {
  if (is.null(callback)) {
    return(NULL)
  }
  if (!is.function(callback)) {
    stop("Expected a function for 'callback'")
  }
  if (length(formals(callback)) != 3) {
    stop("Expected a function with 3 arguments for 'callback'")
  }
  list(callback, new.env(parent = environment(callback)))
}


DOPRI_IDX_TARGET <- 1L
DOPRI_IDX_PARMS <- 2L
DOPRI_IDX_ENV <- 3L
DOPRI_IDX_OUTPUT <- 4L
DOPRI_IDX_EVENT <- 5L

## Must match up with enum dopri_verbose in dopri.h
VERBOSE_QUIET <- 0L
VERBOSE_STEP <- 1L
VERBOSE_EVAL <- 2L
