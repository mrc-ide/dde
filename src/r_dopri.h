#include <R.h>
#include <Rinternals.h>
#include "dopri.h"

// Main interface:
SEXP r_dopri(SEXP r_y_initial, SEXP r_times, SEXP r_func, SEXP r_data,
             SEXP r_n_out, SEXP r_output, SEXP r_data_is_real,
             SEXP r_rtol, SEXP r_atol,
             SEXP r_step_size_min, SEXP r_step_size_max,
             SEXP r_step_size_initial, SEXP r_step_max_n,
             // Critical times and events:
             SEXP r_tcrit, SEXP r_is_event, SEXP r_events,
             // Other:
             SEXP r_use_853,
             SEXP r_stiff_check,
             SEXP r_verbose, SEXP r_callback,
             // Return information
             SEXP r_n_history, SEXP r_grow_history, SEXP r_return_history,
             SEXP r_return_initial, SEXP r_return_statistics,
             SEXP r_return_pointer);
SEXP r_dopri_continue(SEXP r_ptr, SEXP r_y_initial, SEXP r_times,
                      SEXP r_data, SEXP r_data_is_real, SEXP r_tcrit,
                      // Return information:
                      SEXP r_return_history, SEXP r_return_initial,
                      SEXP r_return_statistics, SEXP r_return_pointer);
SEXP r_dopri_copy(SEXP r_ptr);
SEXP r_ylag(SEXP r_t, SEXP r_idx);

void dde_r_harness(size_t n, double t, const double *y, double *dydt,
                   const void *data);
void dde_r_output_harness(size_t n, double t, const double *y,
                          size_t n_out, double *out, const void *data);
void dde_r_event_harness(size_t n, double t, double *y, void *data);

// Internal
SEXP dopri_ptr_create(dopri_data *obj);
void dopri_ptr_finalizer(SEXP extPtr);
dopri_data* dopri_ptr_get(SEXP r_ptr);
void r_dopri_error(dopri_data* obj);
void r_dopri_cleanup(dopri_data *obj, SEXP r_ptr, SEXP r_y,
                     bool return_history, bool return_statistics,
                     bool return_pointer);
