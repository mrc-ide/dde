#include "dopri5.h"
#include <R.h>
#include <Rinternals.h>

// There are more arguments from lsoda not implemented here that will
// be needed:
//
// - rootfunc (once root-finding is supported)
// - verbose (some sort of verbose output will be useful)
// - nroot (not sure)
// - hmin, hmax, hini
// - ynames
// - maxsteps
// - nout
// - outnames
// - forcings
// - events
//
// Some of these are big issues, some are small!
SEXP r_dopri5(SEXP r_y_initial, SEXP r_times, SEXP r_func, SEXP r_data,
              SEXP r_n_out, SEXP r_output,
              SEXP r_rtol, SEXP r_atol, SEXP r_data_is_real,
              SEXP r_tcrit,
              SEXP r_n_history, SEXP r_return_history,
              SEXP r_return_initial, SEXP r_return_statistics) {
  double *y_initial = REAL(r_y_initial);
  size_t n = length(r_y_initial);

  size_t n_times = LENGTH(r_times);
  double *times = REAL(r_times);

  size_t n_tcrit = 0;
  double *tcrit = NULL;
  if (r_tcrit != R_NilValue) {
    n_tcrit = LENGTH(r_tcrit);
    tcrit = REAL(r_tcrit);
  }

  deriv_func *func = (deriv_func*)R_ExternalPtrAddr(r_func);
  void *data = NULL;
  if (TYPEOF(r_data) == REALSXP && INTEGER(r_data_is_real)[0]) {
    data = (void*) REAL(r_data);
  } else if (TYPEOF(r_data) == EXTPTRSXP) {
    data = R_ExternalPtrAddr(r_data);
  } else {
    data = (void*) r_data;
  }

  size_t n_history = (size_t)INTEGER(r_n_history)[0];
  bool return_history = INTEGER(r_return_history)[0];
  bool return_initial = INTEGER(r_return_initial)[0];
  bool return_statistics = INTEGER(r_return_statistics)[0];
  size_t nt = return_initial ? n_times : n_times - 1;

  size_t n_out = INTEGER(r_n_out)[0];
  output_func *output = NULL;
  double *out = NULL;
  SEXP r_out = R_NilValue;
  if (n_out > 0) {
    output = (output_func*)R_ExternalPtrAddr(r_output);
    r_out = PROTECT(allocMatrix(REALSXP, n_out, nt));
    out = REAL(r_out);
  }

  // TODO: as an option save the conditions here.  That's not too bad
  // because we just don't pass through REAL(r_y) but REAL(r_y) +
  // n.  We do have to run the output functions once more though.
  dopri5_data* obj = dopri5_data_alloc(func, n, output, n_out, data, n_history);
  obj->rtol = REAL(r_rtol)[0];
  obj->atol = REAL(r_atol)[0];

  SEXP r_y = PROTECT(allocMatrix(REALSXP, n, nt));

  double *y = REAL(r_y);
  if (return_initial) {
    memcpy(y, y_initial, n * sizeof(double));
    y += n;
    if (n_out > 0) {
      obj->output(obj->n, times[0], y_initial, obj->n_out, out, obj->data);
      out += n_out;
    }
  }

  dopri5_integrate(obj, y_initial, times, n_times, tcrit, n_tcrit, y, out);

  if (obj->error) {
    Rf_error("Integration failure with code: %d", obj->code);
  }

  if (return_history) {
    size_t nh = ring_buffer_used(obj->history, 0);
    SEXP history = PROTECT(allocMatrix(REALSXP, obj->history_len, nh));
    ring_buffer_memcpy_from(REAL(history), obj->history, nh);
    setAttrib(r_y, install("history"), history);
    UNPROTECT(1);
  }

  if (n_out > 0) {
    setAttrib(r_y, install("output"), r_out);
    UNPROTECT(1);
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
    UNPROTECT(2);
  }

  dopri5_data_free(obj);

  UNPROTECT(1);
  return r_y;
}

void dde_r_harness(size_t n, double t, double *y, double *dydt, void *data) {
  SEXP d = (SEXP)data;
  SEXP
    target = VECTOR_ELT(d, 0),
    parms = VECTOR_ELT(d, 1),
    rho = VECTOR_ELT(d, 2);
  SEXP r_t = PROTECT(ScalarReal(t));
  SEXP r_y = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(r_y), y, n * sizeof(double));
  SEXP call = PROTECT(lang4(target, r_t, r_y, parms));
  SEXP ans = PROTECT(eval(call, rho));
  memcpy(dydt, REAL(ans), n * sizeof(double));
  UNPROTECT(4);
}

void dde_r_output_harness(size_t n, double t, double *y,
                          size_t n_out, double *out, void *data) {
  SEXP d = (SEXP)data;
  SEXP
    parms = VECTOR_ELT(d, 1),
    rho = VECTOR_ELT(d, 2),
    output = VECTOR_ELT(d, 3);
  SEXP r_t = PROTECT(ScalarReal(t));
  SEXP r_y = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(r_y), y, n * sizeof(double));
  SEXP call = PROTECT(lang4(output, r_t, r_y, parms));
  SEXP ans = PROTECT(eval(call, rho));
  memcpy(out, REAL(ans), n_out * sizeof(double));
  UNPROTECT(4);
}
