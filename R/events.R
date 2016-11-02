## TODO: I can either take a list here (time, event) or I can take the
## two things separately.  I'm undecided what is nicest!
##
## This function is destined to be internal, but will be used by
## rlsoda so the design matters quite a lot.
##
## Getting the multiple case working here is a bit of a faff; we'll
## need to ensure that the functions have the right type, length,
## convert to a list, tidy away the R function and create a closure
## (if appropriate); it's all a bit much really!  For now, let's just
## do the single case.
check_events <- function(time, event, tcrit = NULL, dllname = "") {
  if (!is.null(tcrit)) {
    assert_numeric(tcrit)
    tcrit <- as.numeric(tcrit)
  }
  if (!is.null(time)) {
    assert_numeric(time)
    time <- as.numeric(time)
  }

  ## Early exit if there are no event:
  if (is.null(event)) {
    if (!is.null(time)) {
      stop("confused") # TODO: proper error
    }
    return(list(event = NULL, is_event = NULL, tcrit = tcrit))
  }

  if (is.null(tcrit)) {
    tcrit <- time
    is_event <- rep(TRUE, length(tcrit))
  } else {
    ## NOTE: We do want to preserve duplicated times within 'time',
    ## but filter out those that interact with tcrit.  The other way
    ## of doing this would be simply to sort all the tcrit times into
    ## time with no event.
    tcrit <- sort(c(setdiff(tcrit, time), time))
    is_event <- tcrit %in% time
  }

  get_event_function <- function(event) {
    event <- find_function_address(event, dllname)

    ## This is the same pattern as for resolving output.  It may be
    ## relaxed in future.
    if (is_r_target) {
      if (!is.function(event)) {
        stop("event must be an R function")
      }
      event <- find_function_address("dde_r_event_harness", "dde")
    } else {
      if (is.function(event)) {
        stop("event must be a compiled function (name or address)")
      }
    }
    event
  }

  is_r_target <- !nzchar(dllname)

  event_function <- if (is_r_target) event else NULL

  ## This should allow passing in:
  ##   function
  ##   NativeSymbolInfo
  ##   pointer
  ##   list of any of the above
  ## and do the right thing
  if (is.list(event) && !inherits(event, "NativeSymbolInfo")) {
    event <- lapply(event, get_event_function)
  } else {
    event <- list(get_event_function(event))
  }

  nt <- length(time)
  ne <- length(event)
  if (ne == nt) {
    if (is_r_target && ne != 1L) {
      event_function <- events_call_and_advance(event_function)
    }
  } else if (ne == 1L) {
    event <- rep(event, nt)
  } else {
    stop("'event' must be a single event or a list of length ", nt)
  }

  list(tcrit = tcrit,
       is_event = is_event,
       event = event,
       event_function = event_function,
       is_r_target = is_r_target)
}


## for the events function
events_call_and_advance <- function(funs) {
  force(funs)
  i <- 0L
  function(...) {
    i <<- i + 1L
    funs[[i]](...)
  }
}
