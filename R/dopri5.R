##' Integrate an ODE or DDE with dopri5.
##'
##' @title Integrate ODE/DDE with dopri5
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
##' @param by_column Logical, indicating if the output should be
##'   returned organised by column (rather than row).  This incurs a
##'   slight cost for transposing the matrices.  If you can work with
##'   matrices that are transposed relative to \code{deSolve}, then
##'   set this to \code{FALSE}.
##'
##' @param ynames Logical, indicating if the output should be named
##'   following the names of the input vector \code{y}.
##'   Alternatively, if \code{ynames} is a character vector of the
##'   same length as \code{y}, these will be used as the output names.
##'
##' @param outnames An optional character vector, used when
##'   \code{n_out} is greater than 0, to name the model output matrix.
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
dopri5 <- function(y, times, func, parms, ...,
                   n_out = 0L, output = NULL,
                   rtol = 1e-6, atol = 1e-6,
                   tcrit = NULL,
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

  assert_scalar(rtol)
  assert_scalar(atol)
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
  } else {
    output <- NULL
  }

  ret <- .Call(Cdopri5, y, times, func, parms,
               n_out, output,
               rtol, atol, parms_are_real, tcrit,
               as.integer(n_history), return_history, return_initial,
               return_statistics)

  if (return_time) {
    ret <- rbind(if (return_initial) times else times[-1],
                 ret, deparse.level = 0)
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

dopri5_interpolate <- function(h, t) {
  nh <- ncol(h)
  nd <- (nrow(h) - 2) / 5
  it <- nrow(h) - 1L
  ih <- nrow(h)

  tr <- c(h[it, 1L], h[it, nh] + h[ih, nh])
  if (min(t) < tr[[1L]] || max(t) > tr[[2L]]) {
    stop("Time falls outside of range of known history [%s, %s]",
         tr[[1L]], tr[[2L]])
  }

  ## Next, we need to put history into the right shape.
  history <- array(h[seq_len(5L * nd), ], c(nd, 5L, nh))
  ## Then pulling apart into the 5 coefficients
  h1 <- history[, 1L, ]
  h2 <- history[, 2L, ]
  h3 <- history[, 3L, ]
  h4 <- history[, 4L, ]
  h5 <- history[, 5L, ]

  idx <- findInterval(t, h[it, ])
  theta  <- (t - h[it, idx]) / h[ih, idx]
  theta1 <- 1.0 - theta

  ret <- matrix(NA_real_, length(t), nd)
  for (i in seq_len(nd)) {
    ret[, i] = h1[i, idx] + theta *
      (h2[i, idx] + theta1 *
       (h3[i, idx] + theta *
        (h4[i, idx] + theta1 *
         h5[i, idx])))
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
