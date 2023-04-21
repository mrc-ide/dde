#ifndef _DIFEQ_H_
#define _DIFEQ_H_

#include <dde/dde.h>
#include <stdbool.h>
#include <stddef.h>
#include <ring/ring.h>
#define R_NO_REMAP 1

// This is *very* similar to the format used by deSolve and dopri, but
// with a couple of differences:
//
//   double time has been replaced by size_t step
//
//   const void *data is passed through (vs deSolve)
//
//   output is always passed through (vs dopri)
//
typedef void difeq_target(size_t n_eq, size_t step,
                          const double *y, double *ynext,
                          size_t n_out, double *output, const void *data);

typedef struct {
  difeq_target * target;
  const void *data;

  size_t n;     // number of equations
  size_t n_out; // number of output variables (possibly zero)

  size_t step0; // initial step number; used in models with delays
  size_t step;  // current step number
  size_t step1; // final step number

  // Steps for integration to report at
  const size_t *steps; // Set of steps to stop at
  size_t n_steps;
  size_t steps_idx;

  double * y0; // initial state
  double * y1; // final state

  // History (if needed)
  //
  // laid out as {idx, t, y, out}
  size_t history_n; // number of elements
  size_t history_len; // element length
  ring_buffer *history;
  size_t history_idx_step;
  size_t history_idx_y;
  size_t history_idx_out;
} difeq_data;

difeq_data* difeq_data_alloc(difeq_target * target,
                             size_t n, size_t n_out, const void *data,
                             size_t n_history, bool grow_history);
void difeq_data_reset(difeq_data *obj, const double *y,
                      const size_t *steps, size_t n_steps);
difeq_data* difeq_data_copy(const difeq_data* obj);
void difeq_data_free(difeq_data *obj);
void difeq_run(difeq_data *obj, const double *y,
               const size_t *steps, size_t n_steps,
               double *y_out, double *out,
               bool return_initial);

size_t get_current_problem_size_difeq(void);

#endif
