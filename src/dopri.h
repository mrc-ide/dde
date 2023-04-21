#ifndef _DOPRI_H_
#define _DOPRI_H_

#include <dde/dde.h>
#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <stddef.h>
#include <ring/ring.h>

// The big thing that is missing for testing this is the ability to
// output anything.  Let's simplify this entirely at the moment and go
// with:
//
// * fixed all dense output
// * no full output history saving (yet)
// * a vector of times to output at, saving into a big array (this is
//   probably the best thing to do).
//
// We'll have to do some extra tricks for general output; it might
// involve actually running the target function an additional time per
// sample to fill in the output variables, but that's not the end of
// the world.
typedef enum return_code {
  ERR_ZERO_TIME_DIFFERENCE =-1,
  ERR_INCONSISTENT_TIME = -2,
  ERR_TOO_MANY_STEPS=-3,
  ERR_STEP_SIZE_TOO_SMALL=-4,
  ERR_STEP_SIZE_VANISHED=-5,
  ERR_YLAG_FAIL=-6,
  ERR_STIFF=-7,
  NOT_SET=0,
  OK_COMPLETE = 1,
  OK_INTERRUPTED = 2 // TODO: never used!
} return_code;

typedef enum dopri_method {
  DOPRI_5,
  DOPRI_853
} dopri_method;

// must match up with the R constants in dopri.R
typedef enum dopri_verbose {
  VERBOSE_QUIET = 0,
  VERBOSE_STEP = 1,
  VERBOSE_EVAL = 2
} dopri_verbose;

typedef void deriv_func(size_t n_eq, double t, const double *y, double *dydt,
                        const void *data);
typedef void output_func(size_t n_eq, double t, const double *y,
                         size_t n_out, double *out, const void *data);
typedef void event_func(size_t n_eq, double t, double *y, void *data);

typedef struct {
  deriv_func* target;  // core rhs function
  output_func* output; // optional output function

  void* data; // model data

  dopri_method method; // switch between (4)5 and 8(5(3))
  size_t order;        // order of integration, based on method

  size_t n;     // number of equations
  size_t n_out; // number of output variables (possibly zero)

  double t0; // initial time; used in models with delays
  double t;  // current time

  dopri_verbose verbose; // be verbose when running
  SEXP callback; // callback function/environment pair (or R_NilValue)

  // Times for integration to report at
  const double *times; // Set of times to stop at
  size_t n_times;
  size_t times_idx;

  // Times for integration to stop at
  const double *tcrit;
  size_t n_tcrit;
  size_t tcrit_idx;

  // Event support
  const bool *is_event;
  event_func **events;

  double * y0; // initial state
  double * y;  // current state
  double * y1; // next state

  // internal state for the dopri step
  double **k;

  // For the delay equation version this is a circular buffer, with
  // each element being length:
  //
  //     5 * n + 2 (DOPRI_5)
  //     8 * n + 2 (DOPRI_853)
  //
  // and with n_history elements.
  //
  // The memory layout of each element;
  //
  // DOPRI5_:
  //
  //   double[n]: y
  //   double[n]: ydiff
  //   double[n]: bspl
  //   double[n]: expr (a function of the above)
  //   double[n]: dens (computed from the ks)
  //   double: t (at 5 * n)
  //   double: h (at 5 * n + 1)
  //
  // DOPRI_853:
  //
  //   double[n]: y
  //   double[n]: ydiff
  //   double[n]: bspl
  //   double[n]: expr (a function of the above)
  //   double[n]: }
  //   double[n]: } four expressions that are written to twice
  //   double[n]: }
  //   double[n]: }
  //   double: t (at 8 * n)
  //   double: h (at 8 * n + 1)
  //
  size_t history_len;
  ring_buffer *history;
  size_t history_idx_time;

  double sign;

  // Control parameters:
  double atol;
  double rtol;
  double step_factor_safe;
  double step_factor_min;
  double step_factor_max;
  double step_size_min;
  double step_size_max;
  double step_size_initial;
  size_t step_max_n; // max number of steps (100000)
  bool step_size_min_allow;
  double step_beta;
  double step_constant; // internal

  // Error reporting
  bool error;
  return_code code;

  // Work statistics
  size_t n_eval;
  size_t n_step;
  size_t n_accept;
  size_t n_reject;

  // Stiffness detection
  size_t stiff_check;
  size_t stiff_n_stiff;
  size_t stiff_n_nonstiff;
} dopri_data;

dopri_data* dopri_data_alloc(deriv_func* target, size_t n,
                             output_func* output, size_t n_out,
                             void *data,
                             dopri_method method,
                             size_t n_history, bool grow_history,
                             dopri_verbose verbose, SEXP callback);
void dopri_data_reset(dopri_data *obj, const double *y,
                      const double *times, size_t n_times,
                      const double *tcrit, size_t n_tcrit,
                      // TODO: naming here may change:
                      const bool *is_event, event_func **events);
dopri_data* dopri_data_copy(const dopri_data* obj);
void dopri_data_free(dopri_data *obj);
void dopri_integrate(dopri_data *obj, const double *y,
                     const double *times, size_t n_times,
                     const double *tcrit, size_t n_tcrit,
                     const bool *is_event, event_func **events,
                     double *y_out, double *out,
                     bool return_initial);

// Wrappers around the two methods:
void dopri_step(dopri_data *obj, double h);
double dopri_error(dopri_data *obj);
void dopri_save_history(dopri_data *obj, double h);
bool dopri_test_stiff(dopri_data *obj, double h);

double dopri_h_new(dopri_data *obj, double fac_old, double h, double err);
double dopri_h_init(dopri_data *obj);

double dopri_interpolate_1(const double *history, dopri_method method,
                           size_t n, double t, size_t i);
void dopri_interpolate_all(const double *history, dopri_method method,
                           size_t n, double t, double *y);
void dopri_interpolate_idx(const double *history, dopri_method method,
                           size_t n, double t, const size_t * idx, size_t nidx,
                           double *y);
void dopri_interpolate_idx_int(const double *history, dopri_method method,
                               size_t n, double t, const int *idx, size_t nidx,
                               double *y);

void dopri_eval(dopri_data *obj, double t, double *y, double *dydt);
void dopri_print_step(dopri_data *obj, double h);
void dopri_print_eval(dopri_data *obj, double t, double *y);
void dopri_callback(dopri_data *obj, double t, double h, double *y);

// Helper
size_t get_current_problem_size_dde(void);
double square(double x);
size_t min_size(size_t a, size_t b);

#endif
