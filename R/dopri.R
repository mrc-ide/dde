##' Integrate an ODE or DDE with dopri.
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
##'   \code{size_t n, double t, double *y, double *dydt, void *data}.
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
##'   throw an error, stopping the integration (note that this differs
##'   from the treatment of \code{hmin} in \code{deSolve::lsoda}).
##'
##' @param step_size_max The largest step size.  By default there is
##'   no maximum step size (Inf) so the solver can take as large a
##'   step as it wants to.  If you have short events you want the
##'   solver to notice, then specify a smaller maximim step size here
##'   (or use \code{tcrit} below).
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
##' @param method The integration method to use, as a string.  The
##'   supported methods are \code{"dopri5"} (5th order method with 4th
##'   order dense output) and \code{"dopri853"} (8th order method with
##'   7th order output and embedded 5th and 3rd order schemes).
##'   Alternatively, use the functions \code{dopri5} or
##'   \code{dopri853} which simply sets this argument.
##'
##' @param n_history Number of history points to retain.  This needs
##'   to be greater than zero for delay differential equations to
##'   work.  Alternatively, this may be greater than zero to return
##'   model outputs that can be inspected later.
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
##' @param by_column Logical, indicating if the output should be
##'   returned organised by column (rather than row).  This incurs a
##'   slight cost for transposing the matrices.  If you can work with
##'   matrices that are transposed relative to \code{deSolve}, then
##'   set this to \code{FALSE}.
##'
##' @param return_initial Logical, indicating if the output should
##'   include the initial conditions (like deSolve).
##'
##' @param return_statistics Logical, indicating if statistics about
##'   the run should be included.  If \code{TRUE}, then an integer
##'   vector containing the number of target evaluations, steps,
##'   accepted steps and rejected steps is returned (the vector is
##'   named).
##'
##' @param return_time Logical, indicating if a row (or column if
##'   \code{by_column} is \code{TRUE}) representing time is included
##'   (this matches deSolve).
##'
##' @param deSolve_compatible Logical, indicating if we should run in
##'   "deSolve compatible" output mode.  This enables the options
##'   \code{by_column}, \code{return_initial} and \code{return_time}.
##'   This affects only some aspects of the output, and not the
##'   calculations themselves.
##'
##' @return At present the return value is transposed relative to
##'   deSolve.  This might change in future.
##'
##' @export
##' @useDynLib dde, .registration = TRUE
dopri <- function(y, times, func, parms, ...,
                  n_out = 0L, output = NULL,
                  rtol = 1e-6, atol = 1e-6,
                  step_size_min = 0, step_size_max = Inf,
                  step_size_initial = 0, step_max_n = 100000L,
                  tcrit = NULL,
                  method = "dopri5",
                  n_history = 0, return_history = n_history > 0, dllname = "",
                  parms_are_real = TRUE,
                  ynames = TRUE, outnames = NULL,
                  by_column = FALSE, return_initial = FALSE,
                  return_statistics=FALSE, return_time = FALSE,
                  deSolve_compatible = FALSE) {
  ## TODO: include "deSolve" mode where we do the transpose, add the
  ## time column too?
  DOTS <- list(...)
  if (length(DOTS) > 0L) {
    stop("Invalid dot arguments!")
  }
  if (deSolve_compatible) {
    by_column <- TRUE
    return_initial <- TRUE
    return_time <- TRUE
  }

  func <- find_function_address(func, dllname)

  is_r_target <- is.function(func)

  if (is_r_target) {
    parms_are_real <- FALSE
    parms <- list(func = func, parms = parms, rho = parent.frame())
    func <- find_function_address("dde_r_harness", "dde")
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
  assert_size(n_history)
  assert_scalar_logical(return_history)
  assert_scalar_logical(parms_are_real)
  assert_scalar_logical(by_column)
  assert_scalar_logical(return_initial)
  assert_scalar_logical(return_statistics)

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

  if (return_time && !is.null(ynames)) {
    ynames <- c("time", ynames)
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

  if (n_out > 0L) {
    output <- find_function_address(output, dllname)
    if (is_r_target) {
      if (!is.function(output)) {
        stop("output must be an R function")
      }
      parms <- c(parms, output)
      output <- find_function_address("dde_r_output_harness", "dde")
    } else {
      if (is.function(output)) {
        stop("output must be a compiled function (name or address)")
      }
    }
    ## Here, if fun is an R function we need to be careful...
  } else if (!is.null(output)) {
    stop("If 'output' is given, then n_out must be specified")
  }

  ret <- .Call(Cdopri, y, as.numeric(times), func, parms,
               as.integer(n_out), output, parms_are_real,
               ## Tolerance:
               rtol, atol,
               ## Step control:
               step_size_min, step_size_max,
               step_size_initial, as.integer(step_max_n),
               ## Other:
               tcrit, use_853,
               ## Return information:
               as.integer(n_history), return_history,
               return_initial, return_statistics)

  if (return_time) {
    at <- attributes(ret)
    ret <- rbind(if (return_initial) times else times[-1],
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

##' @export
##' @rdname dopri
dopri5 <- function(y, times, func, parms, ...) {
  dopri(y, times, func, parms, ..., method="dopri5")
}
##' @export
##' @rdname dopri
dopri853 <- function(y, times, func, parms, ...) {
  dopri(y, times, func, parms, ..., method="dopri853")
}

dopri_interpolate <- function(h, t) {
  nd <- attr(h, "n") # number of equations
  nh <- ncol(h)      # number of history entries
  it <- nrow(h) - 1L # time index
  ih <- nrow(h)      # step size index
  order <- (nrow(h) - 2) / nd # order of integration

  tr <- c(h[it, 1L], h[it, nh] + h[ih, nh])
  if (min(t) < tr[[1L]] || max(t) > tr[[2L]]) {
    stop("Time falls outside of range of known history [%s, %s]",
         tr[[1L]], tr[[2L]])
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
      ret[, i] = h1[i, idx] + theta *
        (h2[i, idx] + theta1 *
         (h3[i, idx] + theta *
          (h4[i, idx] + theta1 *
           h5[i, idx])))
    }
  } else if (order == 8) {
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
  } else {
    stop("Corrupt history object")
  }

  ret
}

find_function_address <- function(fun, dllname = "") {
  if (is.character(fun)) {
    fun <- getNativeSymbolInfo(fun, dllname)$address
  } else if (inherits(fun, "NativeSymbolInfo")) {
    fun <- fun$address
  } else if (!(inherits(fun, "externalptr") || is.function(fun))) {
    stop("Invalid input for 'fun'")
  }
  fun
}

ylag <- function(t, i = NULL) {
  .Call(Cylag, t, i)
}

dopri_methods <- function() {
  c("dopri5", "dopri853")
}