#include <R.h>
#include <Rinternals.h>

SEXP r_dopri(SEXP r_y_initial, SEXP r_times, SEXP r_func, SEXP r_data,
             SEXP r_n_out, SEXP r_output, SEXP r_data_is_real,
             SEXP r_rtol, SEXP r_atol,
             SEXP r_step_size_min, SEXP r_step_size_max,
             SEXP r_step_size_initial, SEXP r_step_max_n,
             SEXP r_tcrit,
             SEXP r_use_853,
             SEXP r_n_history, SEXP r_return_history,
             SEXP r_return_initial, SEXP r_return_statistics,
             SEXP r_return_pointer);
SEXP r_ylag(SEXP r_t, SEXP r_idx);
