#include "dopri5.h"
#include <R.h>
#include <Rinternals.h>

// There are more arguments from lsoda not implemented here that will
// be needed:
//
// - rootfunc (once root-finding is supported)
// - verbose (some sort of verbose output will be useful)
// - nroot (not sure)
// - tcrit (this one is going to be needed)
// - hmin, hmax, hini
// - ynames
// - maxsteps
// - nout
// - outnames
// - forcings
// - events
//
// Some of these are big issues, some are small!
SEXP r_dopri5(SEXP r_y, SEXP r_times, SEXP r_func, SEXP r_data,
              SEXP r_n_out, SEXP r_output,
              SEXP r_rtol, SEXP r_atol, SEXP data_is_real,
              SEXP r_n_history, SEXP r_keep_history) {
  size_t n = length(r_y);
  double *y = REAL(r_y);

  size_t n_times = LENGTH(r_times);
  double *times = REAL(r_times);

  deriv_func *func = (deriv_func*)R_ExternalPtrAddr(r_func);
  void *data = INTEGER(data_is_real)[0] ? (void*) REAL(r_data) : (void*) r_data;
  size_t n_history = (size_t)INTEGER(r_n_history)[0];
  bool keep_history = INTEGER(r_keep_history)[0];

  size_t n_out = INTEGER(r_n_out)[0];
  output_func *output = NULL;
  double *out = NULL;
  SEXP r_out = R_NilValue;
  if (n_out > 0) {
    output = (output_func*)R_ExternalPtrAddr(r_output);
    r_out = PROTECT(allocMatrix(REALSXP, n_out, n_times - 1));
    out = REAL(r_out);
  }

  dopri5_data* obj = dopri5_data_alloc(func, n, output, n_out, data, n_history);
  obj->rtol = REAL(r_rtol)[0];
  obj->atol = REAL(r_atol)[0];

  SEXP ret_y = PROTECT(allocMatrix(REALSXP, n, n_times - 1));
  dopri5_integrate(obj, y, times, n_times, REAL(ret_y), out);

  if (obj->error) {
    Rf_error("Integration failure with code: %d", obj->code);
  }

  if (keep_history) {
    size_t nh = ring_buffer_used(obj->history, 0);
    SEXP history = PROTECT(allocMatrix(REALSXP, obj->history_len, nh));
    ring_buffer_memcpy_from(REAL(history), obj->history, nh);
    setAttrib(ret_y, install("history"), history);
    UNPROTECT(1);
  }

  if (n_out > 0) {
    setAttrib(ret_y, install("output"), r_out);
    UNPROTECT(1);
  }

  UNPROTECT(1);
  return ret_y;
}
