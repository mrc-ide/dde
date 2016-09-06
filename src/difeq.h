#ifndef _DIFEQ_H_
#define _DIFEQ_H_

#include <dde/dde.h>
#include <stdbool.h>
#include <stddef.h>
#include <ring/ring.h>
#define R_NO_REMAP 1

typedef enum return_code {
  ERR_ZERO_STEP_DIFFERENCE = -1,
  ERR_INCONSISTENT_STEP = -2,
  ERR_YLAG_FAIL = -4,
  NOT_SET = 0,
  OK_COMPLETE = 1
} return_code;

// This is *very* similar to the format used by deSolve and dopri, but
// with a couple of differences:
//
//   double time has been augmented by size_t step
//
//   const void *data is passed through
//
typedef void difeq_target(size_t n_eq, size_t step, double time,
                          const double *y, double *ynext,
                          size_t n_out, double *output, const void *data);

typedef struct {
  difeq_target * target;
  void *data;

  size_t n;     // number of equations
  size_t n_out; // number of output variables (possibly zero)

  size_t i0; // initial step number (always zero?)
  size_t i;  // current step number
  size_t i1; // final step number

  double t0;
  double dt;
  double t;

  // Steps for integration to report at
  size_t *steps; // Set of steps to stop at
  size_t n_steps;
  size_t steps_idx;

  double * y0; // initial state

  // History (if needed)
  //
  // laid out as {idx, (t), y, out}
  size_t history_n; // number of elements
  size_t history_len; // element length
  ring_buffer *history;
  size_t history_idx_step;
  size_t history_idx_time;
  size_t history_idx_y;
  size_t history_idx_out;

  // Error reporting
  bool error;
  return_code code;
} difeq_data;

difeq_data* difeq_data_alloc(difeq_target * target,
                             size_t n, size_t n_out, void *data,
                             size_t n_history);
void difeq_data_reset(difeq_data *obj, double *y,
                      size_t *steps, size_t n_steps,
                      double t0, double dt);
void difeq_data_free(difeq_data *obj);
void difeq_run(difeq_data *obj, double *y,
               size_t *steps, size_t n_steps, double t0, double dt,
               double *y_out, double *out,
               bool return_initial);

size_t get_current_problem_size_difeq();

#endif
