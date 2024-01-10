#include "r_dopri.h"
#include "dopri.h"
#include <R.h>
#include <Rinternals.h>
#include "util.h"

enum dopri_idx {
  DOPRI_IDX_TARGET,
  DOPRI_IDX_PARMS,
  DOPRI_IDX_ENV,
  DOPRI_IDX_OUTPUT,
  DOPRI_IDX_EVENT,
};

// There are more arguments from lsoda not implemented here that will
// be needed:
//
// - rootfunc (once root-finding is supported)
// - nroot (not sure)
// - forcings (not sure now)
// - events
//
// Some of these are big issues, some are small!
//

SEXP r_dopri(SEXP r_y_initial, SEXP r_times, SEXP r_func, SEXP r_data,
             SEXP r_n_out, SEXP r_output, SEXP r_data_is_real,
             // Tolerance:
             SEXP r_rtol, SEXP r_atol,
             // Step size control:
             SEXP r_step_size_min, SEXP r_step_size_max,
             SEXP r_step_size_initial, SEXP r_step_max_n,
             SEXP r_step_size_min_allow,
             // Critical times and events:
             SEXP r_tcrit, SEXP r_is_event, SEXP r_events,
             // Other:
             SEXP r_use_853,
             SEXP r_stiff_check,
             SEXP r_verbose, SEXP r_callback,
             // Return information:
             SEXP r_n_history, SEXP r_grow_history, SEXP r_return_history,
             SEXP r_return_initial, SEXP r_return_statistics,
             SEXP r_return_pointer) {
  double *y_initial = REAL(r_y_initial);
  size_t n = length(r_y_initial);

  size_t n_times = LENGTH(r_times);
  double *times = REAL(r_times);

  if (n_times < 2) {
    Rf_error("At least two times must be given");
  }

  size_t n_tcrit = 0;
  double *tcrit = NULL;
  if (r_tcrit != R_NilValue) {
    n_tcrit = LENGTH(r_tcrit);
    tcrit = REAL(r_tcrit);
  }

  // TODO: I should probably overhaul the tcrit to critical_times or
  // something less obscure.
  //
  // OK, the *simplest* way of implementing this is to indicate which
  // of the critical times are going to also be events.  We could do
  // this bit of interleaving in R, in the C/R wrapper, or in the C
  // function itself.
  //
  // The trick is going to be to do the sort, then take the times that
  // are *just about close enough*; we'll want to collapse times that
  // are within eps of each other I think.  However, all of this only
  // needs doing in the case that *both* tevents and tcrit are given.
  // Otherwise it's fairly trivial.
  //
  // I'm increasingly thinking this will be easiest to process in R.
  bool *is_event = (bool*)R_alloc(n_tcrit, sizeof(bool));
  event_func **events = NULL;
  if (r_is_event != R_NilValue) {
    int *tmp = INTEGER(r_is_event);
    for (size_t i = 0; i < n_tcrit; ++i) {
      is_event[i] = tmp[i];
    }
    // Transient storage here; will be reclaimed at the end
    events = (event_func**)R_alloc(n_tcrit, sizeof(event_func*));
    for (size_t i = 0, j = 0; i < n_tcrit; ++i) {
      if (is_event[i]) {
        SEXP el = VECTOR_ELT(r_events, j++);
        if (el == R_NilValue) {
          events[i] = dde_r_event_harness;
        } else {
          events[i] = (event_func*)ptr_fn_get(el);
          if (events[i] == NULL) {
            Rf_error("Was passed null pointer for events[%d]", (int)j);
          }
        }
      }
    }
  } else {
    for (size_t i = 0; i < n_tcrit; ++i) {
      is_event[i] = false;
    }
  }

  // There is probably a nicer way of doing this, but while we have
  // exactly two methods, this is not too bad.
  dopri_method method = INTEGER(r_use_853)[0] ? DOPRI_853 : DOPRI_5;

  deriv_func *func = NULL;
  if (r_func == R_NilValue) {
    func = dde_r_harness;
  } else {
    func = (deriv_func*)ptr_fn_get(r_func);
    if (func == NULL) {
      Rf_error("Was passed null pointer for 'func'");
    }
  }
  void *data = data_pointer(r_data, r_data_is_real);

  size_t n_history = (size_t)INTEGER(r_n_history)[0];
  bool grow_history = INTEGER(r_grow_history)[0];
  dopri_verbose verbose = INTEGER(r_verbose)[0];
  bool return_history = INTEGER(r_return_history)[0];
  bool return_initial = INTEGER(r_return_initial)[0];
  bool return_statistics = INTEGER(r_return_statistics)[0];
  bool return_pointer = INTEGER(r_return_pointer)[0];
  size_t nt = return_initial ? n_times : n_times - 1;

  size_t n_out = INTEGER(r_n_out)[0];
  output_func *output = NULL;
  double *out = NULL;
  SEXP r_out = R_NilValue;

  if (n_out > 0) {
    if (r_output == R_NilValue) {
      output = dde_r_output_harness;
    } else {
      output = (output_func*)ptr_fn_get(r_output);
      if (output == NULL) {
        Rf_error("Was passed null pointer for 'output'");
      }
    }
  }

  // TODO: as an option save the conditions here.  That's not too bad
  // because we just don't pass through REAL(r_y) but REAL(r_y) +
  // n.  We do have to run the output functions once more though.
  dopri_data* obj = dopri_data_alloc(func, n, output, n_out, data,
                                     method, n_history, grow_history,
                                     verbose, r_callback);
  // This is to prevent leaks in case of early exit.  If we don't make
  // it to the end of the function (for any reason, including an error
  // call in a user function, etc) R will clean up for us once it
  // garbage collects ptr.  Because R resets the protection stack on
  // early exit it is guaranteed to get collected at some point.
  SEXP r_ptr = PROTECT(dopri_ptr_create(obj));

  obj->rtol = REAL(r_rtol)[0];
  obj->atol = REAL(r_atol)[0];
  obj->step_size_min = fmax(fabs(REAL(r_step_size_min)[0]), DBL_EPSILON);
  obj->step_size_max = fmin(fabs(REAL(r_step_size_max)[0]), DBL_MAX);
  obj->step_size_initial = REAL(r_step_size_initial)[0];
  obj->step_max_n = INTEGER(r_step_max_n)[0];
  obj->step_size_min_allow = INTEGER(r_step_size_min_allow)[0];

  obj->stiff_check = INTEGER(r_stiff_check)[0];

  SEXP r_y = PROTECT(allocMatrix(REALSXP, n, nt));
  memset(REAL(r_y), 0, n * nt * sizeof(double));

  if (n_out > 0) {
    r_out = PROTECT(allocMatrix(REALSXP, n_out, nt));
    out = REAL(r_out);
    setAttrib(r_y, install("output"), r_out);
    UNPROTECT(1);
  }

  double *y = REAL(r_y);
  dopri_integrate(obj, y_initial, times, n_times, tcrit, n_tcrit,
                  is_event, events, y, out, return_initial);

  r_dopri_cleanup(obj, r_ptr, r_y,
                  return_history, return_statistics, return_pointer);
  UNPROTECT(2);
  return r_y;
}

SEXP r_dopri_copy(SEXP r_ptr) {
  return dopri_ptr_create(dopri_data_copy(dopri_ptr_get(r_ptr)));
}

// Different interface here:
SEXP r_dopri_continue(SEXP r_ptr, SEXP r_y_initial, SEXP r_times,
                      SEXP r_data, SEXP r_data_is_real, SEXP r_tcrit,
                      // Return information:
                      SEXP r_return_history, SEXP r_return_initial,
                      SEXP r_return_statistics, SEXP r_return_pointer) {
  dopri_data* obj = dopri_ptr_get(r_ptr);
  size_t n = obj->n, n_out = obj->n_out;
  double *y_initial;
  if (r_y_initial == R_NilValue) {
    y_initial = obj->y;
  } else {
    if ((size_t) length(r_y_initial) != n) {
      Rf_error("Incorrect size 'y' on integration restart");
    }
    y_initial = REAL(r_y_initial);
  }

  size_t n_times = LENGTH(r_times);
  double *times = REAL(r_times);
  if (n_times < 2) {
    Rf_error("At least two times must be given");
  }
  if (times[0] != obj->t) {
    Rf_error("Incorrect initial time on integration restart");
  }
  if (obj->sign != copysign(1.0, times[n_times - 1] - times[0])) {
    Rf_error("Incorrect sign for the times");
  }

  // Need to freshly set the data pointer because it could have been
  // garbage collected in the meantime.
  obj->data = data_pointer(r_data, r_data_is_real);

  bool return_history = INTEGER(r_return_history)[0];
  bool return_initial = INTEGER(r_return_initial)[0];
  bool return_statistics = INTEGER(r_return_statistics)[0];
  bool return_pointer = INTEGER(r_return_pointer)[0];
  size_t nt = return_initial ? n_times : n_times - 1;

  size_t n_tcrit = 0;
  double *tcrit = NULL;
  if (r_tcrit != R_NilValue) {
    n_tcrit = LENGTH(r_tcrit);
    tcrit = REAL(r_tcrit);
  }

  SEXP r_y = PROTECT(allocMatrix(REALSXP, n, nt));
  double *y = REAL(r_y);
  SEXP r_out = R_NilValue;
  double *out = NULL;
  if (n_out > 0) {
    r_out = PROTECT(allocMatrix(REALSXP, n_out, nt));
    out = REAL(r_out);
    setAttrib(r_y, install("output"), r_out);
    UNPROTECT(1);
  }

  bool *is_event = (bool*)R_alloc(n_tcrit, sizeof(bool));
  for (size_t i = 0; i < n_tcrit; ++i) {
    is_event[i] = false;
  }

  dopri_integrate(obj, y_initial, times, n_times, tcrit, n_tcrit,
                  is_event, NULL, // FIXME
                  y, out, return_initial);
  r_dopri_cleanup(obj, r_ptr, r_y,
                  return_history, return_statistics, return_pointer);
  UNPROTECT(1);
  return r_y;
}

void r_dopri_error(dopri_data* obj) {
  int code = obj->code;
  double t = obj->t;
  switch (code) {
  case ERR_ZERO_TIME_DIFFERENCE:
    Rf_error("Initialisation failure: Beginning and end times are the same");
    break;
  case ERR_INCONSISTENT_TIME:
    Rf_error("Initialisation failure: Times have inconsistent sign");
    break;
  case ERR_TOO_MANY_STEPS:
    Rf_error("Integration failure: too many steps (at t = %2.5f)", t);
    break;
  case ERR_STEP_SIZE_TOO_SMALL:
    Rf_error("Integration failure: step size too small (at t = %2.5f)", t);
    break;
  case ERR_STEP_SIZE_VANISHED:
    Rf_error("Integration failure: step size vanished (at t = %2.5f)", t);
    break;
  case ERR_YLAG_FAIL:
    if (ring_buffer_size(obj->history, false) == 0) {
      Rf_error("Integration failure: can't use ylag in model with no history");
    } else {
      Rf_error("Integration failure: did not find time in history (at t = %2.5f)", t);
    }
    break;
  case ERR_STIFF:
    Rf_error("Integration failure: problem became stiff (at t = %2.5f)", t);
    break;
  default:
    Rf_error("Integration failure: (code %d) [dde bug]", code); // #nocov
    break;
  }
}  // #nocov

SEXP r_ylag(SEXP r_t, SEXP r_idx) {
  size_t n = get_current_problem_size_dde();
  if (n == 0) {
    Rf_error("Can't call this without being in an integration");
  }
  double t = scalar_double(r_t);
  SEXP r_y;
  // TODO: bad things will happen if negative values of ylag are
  // given!
  if (r_idx == R_NilValue) {
    r_y = PROTECT(allocVector(REALSXP, n));
    ylag_all(t, REAL(r_y));
  } else {
    const size_t ni = length(r_idx);
    r_y = PROTECT(allocVector(REALSXP, ni));
    if (ni == 1) {
      REAL(r_y)[0] = ylag_1(t, r_index(r_idx, n));
    } else {
      ylag_vec(t, r_indices(r_idx, n), ni, REAL(r_y));
    }
  }
  UNPROTECT(1);
  return r_y;
}

void dde_r_harness(size_t n, double t, const double *y, double *dydt,
                   const void *data) {
  SEXP d = (SEXP)data;
  SEXP
    target = VECTOR_ELT(d, DOPRI_IDX_TARGET),
    parms = VECTOR_ELT(d, DOPRI_IDX_PARMS),
    env = VECTOR_ELT(d, DOPRI_IDX_ENV);
  SEXP r_t = PROTECT(ScalarReal(t));
  SEXP r_y = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(r_y), y, n * sizeof(double));
  SEXP call = PROTECT(lang4(target, r_t, r_y, parms));
  SEXP ans = PROTECT(eval(call, env));
  memcpy(dydt, REAL(ans), n * sizeof(double));
  UNPROTECT(4);
}

void dde_r_output_harness(size_t n, double t, const double *y,
                          size_t n_out, double *out, const void *data) {
  SEXP d = (SEXP)data;
  SEXP
    parms = VECTOR_ELT(d, DOPRI_IDX_PARMS),
    env = VECTOR_ELT(d, DOPRI_IDX_ENV),
    output = VECTOR_ELT(d, DOPRI_IDX_OUTPUT);
  SEXP r_t = PROTECT(ScalarReal(t));
  SEXP r_y = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(r_y), y, n * sizeof(double));
  SEXP call = PROTECT(lang4(output, r_t, r_y, parms));
  SEXP ans = PROTECT(eval(call, env));
  memcpy(out, REAL(ans), n_out * sizeof(double));
  UNPROTECT(4);
}

void dde_r_event_harness(size_t n, double t, double *y, void *data) {
  SEXP d = (SEXP)data;
  SEXP
    event = VECTOR_ELT(d, DOPRI_IDX_EVENT),
    parms = VECTOR_ELT(d, DOPRI_IDX_PARMS),
    env = VECTOR_ELT(d, DOPRI_IDX_ENV);
  SEXP r_t = PROTECT(ScalarReal(t));
  SEXP r_y = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(r_y), y, n * sizeof(double));
  SEXP call = PROTECT(lang4(event, r_t, r_y, parms));
  SEXP ans = PROTECT(eval(call, env));
  memcpy(y, REAL(ans), n * sizeof(double));
  SEXP parms_new = getAttrib(ans, install("parms"));
  if (parms_new != R_NilValue) {
    SET_VECTOR_ELT(d, DOPRI_IDX_PARMS, parms_new);
  }
  UNPROTECT(4);
}

void dopri_ptr_finalizer(SEXP r_ptr) {
  void *obj = R_ExternalPtrAddr(r_ptr);
  if (obj) {
    dopri_data_free((dopri_data*) obj);
    R_ClearExternalPtr(r_ptr);
  }
}

SEXP dopri_ptr_create(dopri_data *obj) {
  SEXP r_ptr = PROTECT(R_MakeExternalPtr(obj, R_NilValue, R_NilValue));
  R_RegisterCFinalizer(r_ptr, dopri_ptr_finalizer);
  UNPROTECT(1);
  return r_ptr;
}

dopri_data* dopri_ptr_get(SEXP r_ptr) {
  return (dopri_data*) ptr_get(r_ptr);
}

void r_dopri_cleanup(dopri_data *obj, SEXP r_ptr, SEXP r_y,
                     bool return_history, bool return_statistics,
                     bool return_pointer) {
  if (obj->error) {
    r_dopri_error(obj); // will error
  }

  if (return_history) {
    size_t nh = ring_buffer_used(obj->history, 0);
    SEXP history = PROTECT(allocMatrix(REALSXP, obj->history_len, nh));
    ring_buffer_read(obj->history, REAL(history), nh);
    SEXP r_n = PROTECT(ScalarInteger(obj->n));
    setAttrib(history, install("n"), r_n);
    setAttrib(history, R_ClassSymbol, mkString("dopri_history"));
    setAttrib(r_y, install("history"), history);
    UNPROTECT(2);
  }

  if (return_statistics) {
    SEXP stats = PROTECT(allocVector(INTSXP, 4));
    SEXP stats_nms = PROTECT(allocVector(STRSXP, 4));
    INTEGER(stats)[0] = obj->n_eval;
    SET_STRING_ELT(stats_nms, 0, mkChar("n_eval"));
    INTEGER(stats)[1] = obj->n_step;
    SET_STRING_ELT(stats_nms, 1, mkChar("n_step"));
    INTEGER(stats)[2] = obj->n_accept;
    SET_STRING_ELT(stats_nms, 2, mkChar("n_accept"));
    INTEGER(stats)[3] = obj->n_reject;
    SET_STRING_ELT(stats_nms, 3, mkChar("n_reject"));
    setAttrib(stats, R_NamesSymbol, stats_nms);
    setAttrib(r_y, install("statistics"), stats);
    SEXP r_step_size_initial = PROTECT(ScalarReal(obj->step_size_initial));
    setAttrib(r_y, install("step_size"), r_step_size_initial);
    UNPROTECT(3);
  }

  // Deterministically clean up if we can, otherwise we clean up by R
  // running the finaliser for us when it garbage collects ptr above.
  if (return_pointer) {
    // Need to reset this memory
    obj->times = NULL;
    obj->tcrit = NULL;
    setAttrib(r_y, install("ptr"), r_ptr);
  } else {
    dopri_data_free(obj);
    R_ClearExternalPtr(r_ptr);
  }
}
