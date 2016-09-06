#include <R.h>
#include <Rinternals.h>
#include "dopri.h"

// Main interface:
SEXP r_dopri(SEXP r_y_initial, SEXP r_times, SEXP r_func, SEXP r_data,
             SEXP r_n_out, SEXP r_output, SEXP r_data_is_real,
             SEXP r_rtol, SEXP r_atol,
             SEXP r_step_size_min, SEXP r_step_size_max,
             SEXP r_step_size_initial, SEXP r_step_max_n,
             SEXP r_tcrit,
             SEXP r_use_853,
             SEXP r_stiff_check,
             SEXP r_n_history, SEXP r_return_history,
             SEXP r_return_initial, SEXP r_return_statistics,
             SEXP r_return_pointer);
SEXP r_dopri_continue(SEXP r_ptr, SEXP r_y_initial, SEXP r_times,
                      SEXP r_data, SEXP r_data_is_real, SEXP r_tcrit,
                      // Return information:
                      SEXP r_return_history, SEXP r_return_initial,
                      SEXP r_return_statistics, SEXP r_return_pointer);
SEXP r_dopri_copy(SEXP r_ptr);
SEXP r_ylag(SEXP r_t, SEXP r_idx);

// Internal
SEXP dopri_ptr_create(dopri_data *obj);
void dopri_ptr_finalizer(SEXP extPtr);
dopri_data* dopri_ptr_get(SEXP r_ptr);
void r_dopri_error(dopri_data* obj);
void r_dopri_cleanup(dopri_data *obj, SEXP r_ptr, SEXP r_y, SEXP r_out,
                     bool return_history, bool return_statistics,
                     bool return_pointer);
