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
check_events <- function(event_time, event_function, tcrit = NULL,
                         dllname = "") {
  if (!is.null(tcrit)) {
    assert_numeric(tcrit)
    tcrit <- as.numeric(tcrit)
  }
  if (!is.null(event_time)) {
    assert_numeric(event_time)
    event_time <- as.numeric(event_time)
  }

  ## Early exit if there are no event:
  if (is.null(event_function)) {
    if (!is.null(event_time)) {
      stop("'event_time' given without 'event_function'")
    }
    return(list(event = NULL, is_event = NULL, tcrit = tcrit))
  } else if (is.null(event_time)) {
    stop("'event_function' given without 'event_time'")
  }

  if (is.null(tcrit)) {
    tcrit <- event_time
    is_event <- rep(TRUE, length(tcrit))
  } else {
    ## NOTE: We do want to preserve duplicated times within 'event_time',
    ## but filter out those that interact with tcrit.  The other way
    ## of doing this would be simply to sort all the tcrit times into
    ## event_time with no event_function.
    tcrit <- sort(c(setdiff(tcrit, event_time), event_time))
    is_event <- tcrit %in% event_time
  }

  get_event_function <- function(event) {
    event <- find_function_address(event, dllname)

    ## This is the same pattern as for resolving output.  It may be
    ## relaxed in future.
    if (is_r_target) {
      if (!is.function(event)) {
        stop("'event_function' must be an R function")
      }
      event <- NULL
    } else {
      if (is.function(event)) {
        stop("'event_function' must be a compiled function (name or address)")
      }
    }
    event
  }

  is_r_target <- !nzchar(dllname)

  event_function_r <- if (is_r_target) event_function else NULL

  ## This should allow passing in:
  ##   function
  ##   NativeSymbolInfo
  ##   pointer
  ##   list of any of the above
  ## and do the right thing
  if (is.list(event_function) && !is_native_symbol_info(event_function)) {
    event_function <- lapply(event_function, get_event_function)
  } else {
    event_function <- list(get_event_function(event_function))
  }

  nt <- length(event_time) # not tcrit!
  ne <- length(event_function)
  if (ne == nt) {
    if (is_r_target && ne != 1L) {
      event_function_r <- events_call_and_advance(event_function_r)
    }
  } else if (ne == 1L) {
    event_function <- rep(event_function, nt)
  } else {
    stop("'event_function' must be a single event or a list of length ", nt)
  }

  list(tcrit = tcrit,
       is_event = is_event,
       event = event_function,
       event_r = event_function_r,
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
