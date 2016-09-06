#include <R.h>
#include <Rinternals.h>

SEXP r_difeq(SEXP r_y_initial, SEXP r_steps, SEXP r_target, SEXP r_data,
             SEXP r_n_out,
             SEXP r_t0, SEXP r_dt,
             SEXP r_data_is_real,
             SEXP r_n_history, SEXP r_return_history,
             SEXP r_return_initial);
