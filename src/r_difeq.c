#include "r_difeq.h"
#include "difeq.h"

void r_difeq_throw_error(difeq_data *obj);

SEXP r_difeq(SEXP r_y_initial, SEXP r_steps, SEXP r_target, SEXP r_data,
             SEXP r_n_out,
             SEXP r_t0, SEXP r_dt,
             SEXP r_data_is_real,
             SEXP r_n_history, SEXP r_return_history,
             SEXP r_return_initial) {
  double *y_initial = REAL(r_y_initial);
  size_t n = length(r_y_initial);

  size_t n_steps = LENGTH(r_steps);
  size_t *steps = (size_t*) R_alloc(n_steps, sizeof(size_t));
  // TODO: check no r_steps are negative!
  int* tmp = INTEGER(r_steps);
  for (size_t i = 0; i < n_steps; ++i) {
    if (tmp[i] < 0) {
      Rf_error("Step %d is negative (%d)", i, tmp[i]);
    }
    steps[i] = (size_t) tmp[i];
  }

  double t0 = REAL(r_t0)[0], dt = REAL(r_dt)[0];

  // TODO: check for NULL function pointers here to avoid crashes;
  // also test type?
  difeq_target *target = (difeq_target*)R_ExternalPtrAddr(r_target);
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
  size_t nt = return_initial ? n_steps : n_steps - 1;

  size_t n_out = INTEGER(r_n_out)[0];
  double *out = NULL;
  SEXP r_out = R_NilValue;
  if (n_out > 0) {
    r_out = PROTECT(allocMatrix(REALSXP, n_out, nt));
    out = REAL(r_out);
  }

  // TODO: as an option save the conditions here.  That's not too bad
  // because we just don't pass through REAL(r_y) but REAL(r_y) +
  // n.  We do have to run the output functions once more though.
  difeq_data* obj = difeq_data_alloc(target, n, n_out, data, n_history);

  SEXP r_y = PROTECT(allocMatrix(REALSXP, n, nt));
  double *y = REAL(r_y);

  difeq_run(obj, y_initial, steps, n_steps, t0, dt, y, out, return_initial);

  if (obj->error) {
    r_difeq_throw_error(obj);
  }

  if (return_history) {
    size_t nh = ring_buffer_used(obj->history, 0);
    SEXP history = PROTECT(allocMatrix(REALSXP, obj->history_len, nh));
    ring_buffer_take(obj->history, REAL(history), nh);
    setAttrib(history, install("n"), ScalarInteger(obj->n));
    setAttrib(r_y, install("history"), history);
    UNPROTECT(1);
  }

  if (n_out > 0) {
    setAttrib(r_y, install("output"), r_out);
    UNPROTECT(1);
  }

  difeq_data_free(obj);

  UNPROTECT(1);
  return r_y;
}

void r_difeq_throw_error(difeq_data *obj) {
  int code = obj->code;
  double t = obj->t;

  difeq_data_free(obj);
  switch(code) {
  case ERR_ZERO_STEP_DIFFERENCE:
    Rf_error("Initialisation failure: Beginning and end times are the same");
    break;
  case ERR_INCONSISTENT_STEP:
    Rf_error("Initialisation failure: Times not strictly increasing");
    break;
  case ERR_YLAG_FAIL:
    Rf_error("difeq failure: did not find time in history (at t = %2.5f)", t);
    break;
  default:
    Rf_error("difeq failure: (code %d) [dde bug]", code); // #nocov
    break;
  }
}

void difeq_r_harness(size_t n, size_t i, double t,
                     double *y,  double *ynext,
                     size_t n_out, double *output, void *data) {
  SEXP d = (SEXP)data;
  SEXP
    target = VECTOR_ELT(d, 0),
    parms = VECTOR_ELT(d, 1),
    rho = VECTOR_ELT(d, 2);
  SEXP r_i = PROTECT(ScalarReal(i));
  SEXP r_t = PROTECT(ScalarReal(t));
  SEXP r_y = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(r_y), y, n * sizeof(double));
  SEXP call = PROTECT(lang5(target, r_i, r_t, r_y, parms));
  SEXP ans = PROTECT(eval(call, rho));
  memcpy(ynext, REAL(ans), n * sizeof(double));
  if (n_out > 0) {
    double *ans_output = REAL(getAttrib(ans, install("output")));
    memcpy(output, ans_output, n_out * sizeof(double));
  }
  UNPROTECT(5);
}
