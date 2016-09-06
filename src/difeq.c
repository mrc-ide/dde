#include "difeq.h"
#include <R.h>
#include <stddef.h>

// internal functions:
void difeq_store_time(difeq_data *obj);
void fill_na(double *x, size_t n);

difeq_data* difeq_data_alloc(difeq_target* target,
                             size_t n, size_t n_out, void *data,
                             size_t n_history) {
  difeq_data *ret = (difeq_data*) R_Calloc(1, difeq_data);
  ret->target = target;
  ret->data = data;

  ret->n = n;
  ret->n_out = n_out;

  ret->n_steps = 0;
  ret->steps = NULL;

  // State vectors
  ret->y0 = R_Calloc(n, double); // initial

  // TODO: are sizeof(double) and sizeof(size_t) the same - probably
  // not in general.  All we really need is sizeof(size_t) <=
  // sizeof(double) though, which is probably reasonable.
  ret->history_n = n_history;
  if (n_history > 0) {
    if (sizeof(size_t) > sizeof(double)) {
      Rf_error("difeq bug"); // should never trigger
    }
    ret->history_len = 2 + n + n_out;
    ret->history =
      ring_buffer_create(n_history, ret->history_len * sizeof(double));
    ret->history_idx_step = 0;
    ret->history_idx_time = 1;
    ret->history_idx_y = 2;
    ret->history_idx_out = 2 + n;
  } else {
    ret->history = NULL;
    ret->history_len = 0;
    ret->history_idx_step = 0;
    ret->history_idx_time = 0;
    ret->history_idx_y = 0;
    ret->history_idx_out = 0;
  }

  return ret;
}

void difeq_data_reset(difeq_data *obj, double *y,
                      size_t *steps, size_t n_steps,
                      double t0, double dt) {
  obj->error = false;
  obj->code = NOT_SET;
  memcpy(obj->y0, y, obj->n * sizeof(double));

  obj->i0 = steps[0];
  obj->i  = steps[0];
  obj->i1 = steps[n_steps - 1];

  obj->t0 = t0;
  obj->dt = dt;
  obj->t  = t0;

  obj->n_steps = n_steps;
  obj->steps = (size_t*) R_Realloc(obj->steps, n_steps, size_t);
  memcpy(obj->steps, steps, n_steps * sizeof(size_t));
  obj->steps_idx = 1; // skipping the first step!

  // TODO: I don't check that there is at least one time anywhere in
  // *this* routine, but it is checked in r_difeq which calls
  // this.
  if (steps[n_steps - 1] == steps[0]) {
    obj->error = true;
    obj->code = ERR_ZERO_STEP_DIFFERENCE;
    return;
  }

  for (size_t i = 0; i < n_steps - 1; ++i) {
    if (steps[i + 1] <= steps[i]) { // NOTE: Disallows ties
      obj->error = true;
      obj->code = ERR_INCONSISTENT_STEP;
      return;
    }
  }
}

void difeq_data_free(difeq_data *obj) {
  R_Free(obj->y0);

  if (obj->history) {
    ring_buffer_destroy(obj->history);
  }

  R_Free(obj->steps);
  R_Free(obj);
}

// This is super ugly; see dde/src/dopri.c (variable dde_global_obj)
// for details.
static difeq_data *difeq_global_obj;

// Used to query the problem size safely from the interface.c file
size_t get_current_problem_size_difeq() {
  return difeq_global_obj == NULL ? 0 : difeq_global_obj->n;
}

void difeq_run(difeq_data *obj, double *y,
               size_t *steps, size_t n_steps, double t0, double dt,
               double *y_out, double *out,
               bool return_initial) {
  difeq_data_reset(obj, y, steps, n_steps, t0, dt);
  if (obj->error) {
    return;
  }

  double *y_next, *out_next;

  bool store_y_in_history = obj->history_len > 0;

  // Possibly only set this if the number of history variables is
  // nonzero?  Needs to be set before any calls to target() though.
  if (store_y_in_history) {
    difeq_global_obj = obj;
    difeq_store_time(obj);
    double *h = (double*)obj->history->head;
    memcpy(h + obj->history_idx_y, y, obj->n * sizeof(double));
    fill_na(h + obj->history_idx_out, obj->n_out);
    h = ring_buffer_head_advance(obj->history);
    y_next   = h + obj->history_idx_y;
    out_next = h + obj->history_idx_out;
  }

  // If requested, copy initial conditions into the output space
  if (return_initial) {
    memcpy(y_out, y, obj->n * sizeof(double));
    fill_na(out, obj->n_out);
    y_out += obj->n;
    out += obj->n_out;
  }

  if (!store_y_in_history) {
    y_next = y_out;
    out_next = out;
  }

  while (true) {
    // Try to be clever about where we store things to avoid copies
    // that will always happen.
    //
    // If we are using the ring buffer then storing 'y' anywhere but
    // the ring buffer will suffer a copy, so we might as well store
    // it there.
    //
    // Otherwise we'll use the final output array as a place to store
    // things.
    obj->target(obj->n, obj->i, obj->t, y, y_next, obj->n_out, out_next,
                obj->data);
    obj->i++;
    obj->t += obj->dt;
    y = y_next;

    if (obj->i == obj->steps[obj->steps_idx]) {
      if (store_y_in_history) {
        memcpy(y_out, y_next,   obj->n     * sizeof(double));
        memcpy(out,   out_next, obj->n_out * sizeof(double));
        y_next   = ((double*) obj->history->head) + obj->history_idx_y;
        out_next = ((double*) obj->history->head) + obj->history_idx_out;
      } else {
        y_next += obj->n;
        out_next += obj->n_out;
      }
      y_out += obj->n;
      out += obj->n_out;
      obj->steps_idx++;
    }

    if (store_y_in_history) {
      difeq_store_time(obj);
      double* h = (double*) ring_buffer_head_advance(obj->history);
      y_next = h + obj->history_idx_y;
      out_next = h + obj->history_idx_out;
    }

    if (obj->i == obj->i1) {
      break;
    }
  }

  difeq_global_obj = NULL;
}

// These bits are all nice and don't use any globals
const double* difeq_find_step(difeq_data *obj, size_t step) {
  int offset = obj->i < step;
  const void *h = NULL;
  if (obj->history != NULL && offset >= 0) {
    h = ring_buffer_head_offset(obj->history, (size_t) offset);
  }
  if (h == NULL) {
    obj->error = true;
    obj->code = ERR_YLAG_FAIL;
  }
  return (double*) h;
}

double yprev_1(size_t step, size_t i) {
  if (step <= difeq_global_obj->i0) {
    return difeq_global_obj->y0[i];
  } else {
    const double * h = difeq_find_step(difeq_global_obj, step);
    if (h == NULL) {
      return NA_REAL;
    } else {
      return h[i];
    }
  }
}

void yprev_all(size_t step, double *y) {
  if (step <= difeq_global_obj->i0) {
    memcpy(y, difeq_global_obj->y0, difeq_global_obj->n * sizeof(double));
  } else {
    const double * h = difeq_find_step(difeq_global_obj, step);
    if (h != NULL) {
      memcpy(y, h, difeq_global_obj->n * sizeof(double));
    }
  }
}

void yprev_vec(size_t step, const size_t *idx, size_t nidx, double *y) {
  if (step <= difeq_global_obj->i0) {
    for (size_t i = 0; i < nidx; ++i) {
      y[i] = difeq_global_obj->y0[idx[i]];
    }
  } else {
    const double * h = difeq_find_step(difeq_global_obj, step);
    if (h != NULL) {
      for (size_t i = 0; i < nidx; ++i) {
        y[i] = h[idx[i]];
      }
    }
  }
}

void yprev_vec_int(size_t step, const int *idx, size_t nidx, double *y) {
  if (step <= difeq_global_obj->i0) {
    for (size_t i = 0; i < nidx; ++i) {
      y[i] = difeq_global_obj->y0[idx[i]];
    }
  } else {
    const double * h = difeq_find_step(difeq_global_obj, step);
    if (h != NULL) {
      for (size_t i = 0; i < nidx; ++i) {
        y[i] = h[idx[i]];
      }
    }
  }
}

// Utility
void difeq_store_time(difeq_data *obj) {
  double *h = (double*) obj->history->head;
  h[obj->history_idx_step] = obj->i;
  h[obj->history_idx_time] = obj->t;
}

void fill_na(double *x, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    x[i] = NA_REAL;
  }
}