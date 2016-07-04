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

  assert_scalar(rtol)
  assert_scalar(atol)
  assert_size(n_history)
  assert_scalar_logical(keep_history)
  assert_scalar_logical(parms_are_real)

  assert_size(n_out)
  if (n_out > 0L) {
    output <- find_function_address(output, dllname)
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
  } else if (!inherits(fun, "externalptr")) {
    stop("Invalid input for 'fun'")
  }
  fun
}
