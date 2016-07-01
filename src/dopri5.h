#include <R.h> // dragging in a big include, strip down later
#include <stdbool.h>
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
  ERR_INCONSISTENT=-1,
  ERR_TOO_MANY_STEPS=-2,
  ERR_STEP_SIZE_TOO_SMALL=-3,
  ERR_STIFF=-4,
  ERR_YLAG_ERROR=-5,
  ERR_NOT_INITIALISED=-6,
  NOT_SET=0,
  OK_COMPLETE = 1,
  OK_INTERRUPTED = 2
} return_code;

typedef void deriv_func(int n_eq, double t, double *y, double *dydt,
                        void *data);

typedef struct {
  deriv_func* target;
  void* data;

  size_t n;  // number of equations
  bool initialised;

  double t0; // initial time (needed?)
  double t;  // current time
  double *times; // Set of times to stop at
  size_t n_times;
  size_t times_idx;

  double * y0; // initial state (needed?)
  double * y;  // current state
  double * y1; // next state

  // TODO: If I generalise this over orders, then this becomes a **
  // and the length is stored here -- not too bad.
  double * k1; // internal state for the dopri step
  double * k2;
  double * k3;
  double * k4;
  double * k5;
  double * k6;
  double * ysti;

  // For the dde version this is going to become a circular buffer,
  // but we need a place to keep this here now.  This has length:
  //
  //     5 * n + 2
  //
  // long, (and for dde there will be mxst of these).
  //
  // The memory layout of each element:
  //   double[n]: y
  //   double[n]: deltay
  //   double[n]: bspl
  //   double[n]: expr (a function of the above)
  //   double[n]: dens (computed from the ks)
  //   double: t
  //   double: h
  //
  // Because 'n' is dynamically sized we'll not be able to put these
  // in a circular buffer easily (unless boost's circular buffer is
  // very clever) so I'm going to have do pointer arithmetic here.
  size_t history_len;
  ring_buffer *history;
  size_t history_time_idx;

  double sign;

  // Control parameters:
  double atol;
  double rtol;
  double step_factor_safe;
  double step_factor_min;
  double step_factor_max;
  double step_size_max;
  double step_size_initial;
  size_t step_max_n; // max number of steps (100000)
  double step_beta;

  // Error reporting
  bool error;
  return_code code;

  // Work statistics
  size_t n_eval;
  size_t n_step;
  size_t n_accept;
  size_t n_reject;
} dopri5_data;

dopri5_data* dopri5_data_alloc(deriv_func* target, size_t n, void *data,
                               size_t n_history);
void dopri5_data_reset(dopri5_data *obj, double *y,
                       double *times, size_t n_times);
void dopri5_data_free(dopri5_data *obj);
void dopri5_integrate(dopri5_data *obj, double *y,
                      double *times, size_t n_times,
                      double *y_out);
void dopri5_step(dopri5_data *obj, double h);
double dopri5_error(dopri5_data *obj);
double dopri5_h_new(dopri5_data *obj, double fac_old, double h, double err);
double dopri5_h_init(dopri5_data *obj);

double dopri5_interpolate_1(double *history, size_t n, double t, size_t i);
void dopri5_interpolate_all(double *history, size_t n, double t,
                            double *y, size_t y_stride);
void dopri5_interpolate_idx(double *history, size_t n, double t,
                            size_t * idx, size_t nidx,
                            double *y);
