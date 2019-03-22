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
##' @param ... Additional arguments passed through to
##'   \code{\link{difeq}}.
##'
##' @param as_array Logical, indicating if the output should be
##'   converted into an array.  If \code{TRUE} then \code{res[, , i]}
##'   will contain the \code{i}'th replicate, if \code{FALSE} then
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
difeq_replicate <- function(n, y, ..., as_array = TRUE) {
  ## It's possible that this could be made slightly more efficient by
  ## shifting the allocations all the way into C but that can be done
  ## later preserving this inteface if that seems necessary.

  ## It's not possible to replicate zero times.  It could be done by
  ## predicting the size of the output and returning an appropriately
  ## sized and named array (with zero in the last dimension) but it
  ## seems weird to me.  For the non-simplified case it's not so bad
  ## because we can just return a zero-length list.
  assert_positive_integer(n)

  varying_initial_state <- is.list(y)

  if (varying_initial_state) {
    n_names <- names(y)
    assert_length(y, n)
    if (!all(lengths(y) == length(y[[1L]]))) {
      stop("All 'y' lengths must be the same")
    }
  } else {
    n_names <- NULL
  }

  if (varying_initial_state) {
    res <- lapply(seq_len(n), function(i) difeq(y[[i]], ...))
  } else {
    res <- lapply(seq_len(n), function(i) difeq(y, ...))
  }

  if (as_array) {
    res_array <- array(unlist(res), c(dim(res[[1L]]), n))
    nms <- dimnames(res[[1L]])
    if (!is.null(nms)) {
      dimnames(res_array) <- c(nms, list(n_names))
    }

    ## In theory we could inspect '...' for 'restartable' (or pass
    ## things more explicitly) but that seems harder work given we
    ## exclude n = 0 above.
    res_array <- collect_attributes(res, "restart_data", res_array)
    res_array <- collect_attributes(res, "ptr", res_array)
    res <- res_array
  } else {
    names(res) <- n_names
  }

  res
}
