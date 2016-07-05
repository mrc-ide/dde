##' Integrate an ODE or DDE with dopri5.
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
##' @param n_history Number of history points to retain.  This needs
##'   to be greater than zero for delay differential equations to
##'   work.  Alternatively, this may be greater than zero to return
##'   model outputs that can be inspected later.
##'
##' @param keep_history Logical indicating if history should be
##'   retained or discarded.  By default, history is retained
##'
##' @param dllname Name of the shared library (without extension) to
##'   find the function \code{func} (and \code{output} if given) in
##'   the case where \code{func} refers to compiled function.
##'
##' @param parms_are_real Logical, indicating if \code{parms} should
##'   be treated as vector of doubles by \code{func} (when it is a
##'   compiled function).  If \code{FALSE}, then \code{parms} is
##'   passed through unmodified as a \code{SEXP}.  If \code{TRUE},
##'   then \code{REAL(parms)}, which is \code{double*} is passed
##'   through.
##'
##' @return At present the return value is transposed relative to
##'   deSolve.  This might change in future.
##'
##' @export
##' @useDynLib dde
dopri5 <- function(y, times, func, parms, ...,
                   n_out = 0L, output = NULL,
                   rtol = 1e-6, atol = 1e-6,
                   n_history = 0, keep_history = n_history > 0, dllname = "",
                   parms_are_real = TRUE) {
  DOTS <- list(...)
  if (length(DOTS) > 0L) {
    stop("Invalid dot arguments!")
  }
  ## TODO: will need to support R functions here at some point, for
  ## completeness sake.  This is not that bad to do as they just
  ## become callbacks bound to an environment...
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
  assert_scalar_logical(keep_history)
  assert_scalar_logical(parms_are_real)

  assert_size(n_out)
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

  .Call("r_dopri5", y, times, func, parms,
        n_out, output,
        rtol, atol, parms_are_real,
        as.integer(n_history), keep_history,
        PACKAGE="dde")
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
  .Call("r_ylag", t, i, PACKAGE = "dde")
}
