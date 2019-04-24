#include <R.h>
#include <Rinternals.h>
#include "difeq.h"

// Main interface
SEXP r_difeq(SEXP r_y_initial, SEXP r_steps, SEXP r_target, SEXP r_data,
             SEXP r_n_out,
             SEXP r_data_is_real,
             SEXP r_n_history, SEXP r_grow_history, SEXP r_return_history,
             SEXP r_return_initial, SEXP r_return_pointer);
SEXP r_difeq_continue(SEXP r_ptr, SEXP r_y_initial, SEXP r_steps,
                      SEXP r_data, SEXP r_data_is_real,
                      SEXP r_return_history, SEXP r_return_initial,
                      SEXP r_return_pointer);
SEXP r_difeq_copy(SEXP r_ptr);
SEXP r_yprev(SEXP r_t, SEXP r_idx);

void difeq_r_harness(size_t n, size_t step,
                     const double *y, double *ynext,
                     size_t n_out, double *output, const void *data);

// Internal
SEXP difeq_ptr_create(difeq_data *obj);
void difeq_ptr_finalizer(SEXP extPtr);
difeq_data* difeq_ptr_get(SEXP r_ptr);
void r_difeq_cleanup(difeq_data *obj, SEXP r_ptr, SEXP r_y,
                     bool return_history, bool return_pointer);
