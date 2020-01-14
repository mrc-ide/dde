#include "r_difeq.h"
#include "difeq.h"
#include "util.h"

SEXP r_difeq(SEXP r_y_initial, SEXP r_steps, SEXP r_target, SEXP r_data,
             SEXP r_n_out,
             SEXP r_data_is_real,
             SEXP r_n_history, SEXP r_grow_history, SEXP r_return_history,
             SEXP r_return_initial, SEXP r_return_pointer) {
  double *y_initial = REAL(r_y_initial);
  size_t n = length(r_y_initial);

  size_t n_steps = LENGTH(r_steps);
  // This is required to avoid the issue that int* and size_t* need
  // not be the same.  I should probably just relax this within run
  // instead (TODO)
  size_t *steps = (size_t*) R_alloc(n_steps, sizeof(size_t));
  int* tmp = INTEGER(r_steps);
  for (size_t i = 0; i < n_steps; ++i) {
    steps[i] = (size_t) tmp[i];
  }

  difeq_target *target = NULL;
  if (r_target == R_NilValue) {
    target = difeq_r_harness;
  } else {
    target = (difeq_target*)ptr_fn_get(r_target);
    if (target == NULL) {
      Rf_error("Was passed null pointer for 'target'");
    }
  }
  void *data = data_pointer(r_data, r_data_is_real);

  size_t n_history = (size_t)INTEGER(r_n_history)[0];
  bool return_history = INTEGER(r_return_history)[0];
  bool return_initial = INTEGER(r_return_initial)[0];
  bool return_pointer = INTEGER(r_return_pointer)[0];
  bool grow_history = INTEGER(r_grow_history)[0];
  size_t nt = return_initial ? n_steps : n_steps - 1;

  size_t n_out = INTEGER(r_n_out)[0];
  double *out = NULL;
  SEXP r_out = R_NilValue;

  // TODO: as an option save the conditions here.  That's not too bad
  // because we just don't pass through REAL(r_y) but REAL(r_y) +
  // n.  We do have to run the output functions once more though.
  difeq_data* obj = difeq_data_alloc(target, n, n_out, data,
                                     n_history, grow_history);

  // This is to prevent leaks in case of early exit.  If we don't make
  // it to the end of the function (for any reason, including an error
  // call in a user function, etc) R will clean up for us once it
  // garbage collects ptr.  Because R resets the protection stack on
  // early exit it is guaranteed to get collected at some point.
  SEXP r_ptr = PROTECT(difeq_ptr_create(obj));

  // Then solve things:
  SEXP r_y = PROTECT(allocMatrix(REALSXP, n, nt));
  double *y = REAL(r_y);

  if (n_out > 0) {
    r_out = PROTECT(allocMatrix(REALSXP, n_out, nt));
    out = REAL(r_out);
    setAttrib(r_y, install("output"), r_out);
    UNPROTECT(1);
  }

  GetRNGstate();
  difeq_run(obj, y_initial, steps, n_steps, y, out, return_initial);
  PutRNGstate();

  r_difeq_cleanup(obj, r_ptr, r_y, return_history, return_pointer);

  UNPROTECT(2);
  return r_y;
}

SEXP r_difeq_copy(SEXP r_ptr) {
  return difeq_ptr_create(difeq_data_copy(difeq_ptr_get(r_ptr)));
}

// Different interface here:
SEXP r_difeq_continue(SEXP r_ptr, SEXP r_y_initial, SEXP r_steps,
                      SEXP r_data, SEXP r_data_is_real,
                      // Return information:
                      SEXP r_return_history, SEXP r_return_initial,
                      SEXP r_return_pointer) {
  difeq_data* obj = difeq_ptr_get(r_ptr);
  size_t n = obj->n, n_out = obj->n_out;
  double *y_initial;
  if (r_y_initial == R_NilValue) {
    y_initial = obj->y1;
  } else {
    if ((size_t) length(r_y_initial) != n) {
      Rf_error("Incorrect size 'y' on simulation restart");
    }
    y_initial = REAL(r_y_initial);
  }

  size_t n_steps = LENGTH(r_steps);
  size_t *steps = (size_t*) R_alloc(n_steps, sizeof(size_t));
  int* tmp = INTEGER(r_steps);
  for (size_t i = 0; i < n_steps; ++i) {
    steps[i] = (size_t) tmp[i];
  }
  if (n_steps < 2) {
    Rf_error("At least two steps must be given");
  }
  if (steps[0] != obj->step) {
    Rf_error("Incorrect initial step on simulation restart");
  }

  // Need to freshly set the data pointer because it could have been
  // garbage collected in the meantime.
  obj->data = data_pointer(r_data, r_data_is_real);

  bool return_history = INTEGER(r_return_history)[0];
  bool return_initial = INTEGER(r_return_initial)[0];
  bool return_pointer = INTEGER(r_return_pointer)[0];
  size_t ns = return_initial ? n_steps : n_steps - 1;

  SEXP r_y = PROTECT(allocMatrix(REALSXP, n, ns));
  double *y = REAL(r_y);
  SEXP r_out = R_NilValue;
  double *out = NULL;
  if (n_out > 0) {
    r_out = PROTECT(allocMatrix(REALSXP, n_out, ns));
    out = REAL(r_out);
    setAttrib(r_y, install("output"), r_out);
    UNPROTECT(1);
  }

  GetRNGstate();
  difeq_run(obj, y_initial, steps, n_steps, y, out, return_initial);
  PutRNGstate();

  r_difeq_cleanup(obj, r_ptr, r_y, return_history, return_pointer);

  UNPROTECT(1);
  return r_y;
}

SEXP r_yprev(SEXP r_i, SEXP r_idx) {
  size_t n = get_current_problem_size_difeq();
  if (n == 0) {
    Rf_error("Can't call this without being in an integration");
  }
  int step = scalar_int(r_i);
  SEXP r_y;
  if (r_idx == R_NilValue) {
    r_y = PROTECT(allocVector(REALSXP, n));
    yprev_all(step, REAL(r_y));
  } else {
    const size_t ni = length(r_idx);
    r_y = PROTECT(allocVector(REALSXP, ni));
    if (ni == 1) {
      REAL(r_y)[0] = yprev_1(step, r_index(r_idx, n));
    } else {
      yprev_vec(step, r_indices(r_idx, n), ni, REAL(r_y));
    }
  }
  UNPROTECT(1);
  return r_y;
}

void difeq_r_harness(size_t n, size_t step,
                     const double *y, double *ynext,
                     size_t n_out, double *output, const void *data) {
  SEXP d = (SEXP)data;
  SEXP
    target = VECTOR_ELT(d, 0),
    parms = VECTOR_ELT(d, 1),
    rho = VECTOR_ELT(d, 2);
  SEXP r_step = PROTECT(ScalarInteger(step));
  SEXP r_y = PROTECT(allocVector(REALSXP, n));
  memcpy(REAL(r_y), y, n * sizeof(double));
  SEXP call = PROTECT(lang4(target, r_step, r_y, parms));
  SEXP ans = PROTECT(eval(call, rho));
  // Ensure that we get sensible output from the target function:
  if ((size_t)length(ans) != n) {
    Rf_error("Incorrect length variable output: expected %d, recieved %d",
             n, length(ans));
  }
  memcpy(ynext, REAL(ans), n * sizeof(double));
  if (n_out > 0) {
    SEXP r_output = getAttrib(ans, install("output"));
    if (r_output == R_NilValue) {
      Rf_error("Missing output");
    } else if ((size_t)length(r_output) != n_out) {
      Rf_error("Incorrect length output: expected %d, recieved %d",
               n_out, length(r_output));
    } else if (TYPEOF(r_output) != REALSXP) {
      Rf_error("Incorrect type output");
    }
    memcpy(output, REAL(r_output), n_out * sizeof(double));
  }
  UNPROTECT(4);
}

void difeq_ptr_finalizer(SEXP r_ptr) {
  void *obj = R_ExternalPtrAddr(r_ptr);
  if (obj) {
    difeq_data_free((difeq_data*) obj);
    R_ClearExternalPtr(r_ptr);
  }
}

SEXP difeq_ptr_create(difeq_data *obj) {
  SEXP r_ptr = PROTECT(R_MakeExternalPtr(obj, R_NilValue, R_NilValue));
  R_RegisterCFinalizer(r_ptr, difeq_ptr_finalizer);
  UNPROTECT(1);
  return r_ptr;
}

difeq_data* difeq_ptr_get(SEXP r_ptr) {
  return (difeq_data*) ptr_get(r_ptr);
}

void r_difeq_cleanup(difeq_data *obj, SEXP r_ptr, SEXP r_y,
                     bool return_history, bool return_pointer) {
  if (return_history) {
    size_t nh = ring_buffer_used(obj->history, 0);
    SEXP history = PROTECT(allocMatrix(REALSXP, obj->history_len, nh));
    ring_buffer_read(obj->history, REAL(history), nh);
    SEXP r_n = PROTECT(ScalarInteger(obj->n));
    setAttrib(history, install("n"), r_n);
    setAttrib(r_y, install("history"), history);
    UNPROTECT(2);
  }

  // Deterministically clean up if we can, otherwise we clean up by R
  // running the finaliser for us when it garbage collects ptr above.
  if (return_pointer) {
    obj->steps = NULL;
    setAttrib(r_y, install("ptr"), r_ptr);
  } else {
    difeq_data_free(obj);
    R_ClearExternalPtr(r_ptr);
  }
}
