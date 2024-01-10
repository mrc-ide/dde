#include "difeq.h"
#include <R.h>
#include <stddef.h>

// internal functions:
void difeq_store_time(difeq_data *obj);
void fill_na(double *x, size_t n);

difeq_data* difeq_data_alloc(difeq_target* target,
                             size_t n, size_t n_out, const void *data,
                             size_t n_history, bool grow_history) {
  difeq_data *ret = (difeq_data*) R_Calloc(1, difeq_data);
  overflow_action on_overflow =
    grow_history ? OVERFLOW_GROW : OVERFLOW_OVERWRITE;

  ret->target = target;
  ret->data = data;

  ret->n = n;
  ret->n_out = n_out;

  ret->n_steps = 0;
  ret->steps = NULL;

  // State vectors
  ret->y0 = R_Calloc(n, double); // initial
  ret->y1 = R_Calloc(n, double); // final (only used on restart?)

  ret->history_n = n_history;
  if (n_history > 0) {
    ret->history_len = 1 + n + n_out;
    ret->history =
      ring_buffer_create(n_history, ret->history_len * sizeof(double),
                         on_overflow);
    ret->history_idx_step = 0;
    ret->history_idx_y = 1;
    ret->history_idx_out = 1 + n;
  } else {
    ret->history = NULL;
    ret->history_len = 0;
    ret->history_idx_step = 0;
    ret->history_idx_y = 0;
    ret->history_idx_out = 0;
  }

  return ret;
}

difeq_data* difeq_data_copy(const difeq_data* obj) {
  size_t n_history =
    obj->history == NULL ? 0 : ring_buffer_size(obj->history, false);
  bool grow_history =
    obj->history && obj->history->on_overflow == OVERFLOW_GROW;
  difeq_data* ret = difeq_data_alloc(obj->target, obj->n, obj->n_out,
                                     obj->data, n_history,
                                     grow_history);
  // Then update a few things
  ret->step0 = obj->step0;
  ret->step  = obj->step;
  ring_buffer_mirror(obj->history, ret->history);
  ret->history_idx_step = obj->history_idx_step;

  // NOTE: I'm not copying the steps over because the first thing that
  // happens is always a reset and that deals with this.
  ret->steps = NULL;

  // Copy all of the internal state:
  memcpy(ret->y0, obj->y0, obj->n * sizeof(double));
  memcpy(ret->y1, obj->y1, obj->n * sizeof(double));

  return ret;
}

void difeq_data_reset(difeq_data *obj, const double *y,
                      const size_t *steps, size_t n_steps) {
  memcpy(obj->y0, y, obj->n * sizeof(double));

  obj->n_steps = n_steps;
  obj->steps = steps;
  obj->steps_idx = 1; // skipping the first step!

  // TODO: I don't check that there is at least one time anywhere in
  // *this* routine, but it is checked in r_difeq which calls
  // this.
  if (steps[n_steps - 1] == steps[0]) {
    Rf_error("Initialisation failure: Beginning and end steps are the same");
  }

  for (size_t i = 0; i < n_steps - 1; ++i) {
    if (steps[i + 1] <= steps[i]) { // NOTE: Disallows ties
      Rf_error("Initialisation failure: Steps not strictly increasing");
    }
  }

  // TODO: Possibly need to do some work here.
  if (obj->history == NULL || ring_buffer_is_empty(obj->history)) {
    obj->step0 = steps[0];
  } else {
    // This is the restart condition
    double * h = (double*) ring_buffer_tail(obj->history);
    obj->step0 = h[obj->history_idx_step];
  }
  obj->step = steps[0];
  obj->step1 = steps[n_steps - 1];
}

void difeq_data_free(difeq_data *obj) {
  R_Free(obj->y0);
  R_Free(obj->y1);

  if (obj->history) {
    ring_buffer_destroy(obj->history);
  }

  R_Free(obj);
}

// This is super ugly; see dde/src/dopri.c (variable dde_global_obj)
// for details.
static difeq_data *difeq_global_obj;

// Used to query the problem size safely from the interface.c file
size_t get_current_problem_size_difeq(void) {
  return difeq_global_obj == NULL ? 0 : difeq_global_obj->n;
}

void difeq_run(difeq_data *obj, const double *y,
               const size_t *steps, size_t n_steps,
               double *y_out, double *out,
               bool return_initial) {
  difeq_data_reset(obj, y, steps, n_steps);

  double *y_next = NULL, *out_next = NULL;
  double *write_y = y_out, *write_out = out;

  bool has_output = obj->n_out > 0;
  bool store_next_output = false;
  bool store_y_in_history = obj->history_len > 0;

  if (store_y_in_history) {
    difeq_global_obj = obj;
    bool first_entry = ring_buffer_is_empty(obj->history);
    double *h = (double*) obj->history->head;
    memcpy(h + obj->history_idx_y, y, obj->n * sizeof(double));
    // Not 100% sure if we should do this, but it's probably not that unsafe.
    fill_na(h + obj->history_idx_out, obj->n_out);

    if (first_entry) {
      // No need to store the time as that should be checked ahead of
      // time
      difeq_store_time(obj);
      // Don't think that we should error data out here.
      // And we're off:
      h = ring_buffer_head_advance(obj->history);
    }

    y_next = h + obj->history_idx_y;
    out_next = y_next + obj->n;
  }

  // This is used to detect when the buffer has grown - because this
  // triggers a reallocation we can end up with pointers that are
  // invalidated and we need to update that.  This is only ever
  // accessed when 'store_y_in_history' is true.  Set this *after*
  // running advance above as otherwise it could be immediately
  // invalidated.
  const data_t * buffer_data = store_y_in_history ? obj->history->data : NULL;

  // If requested, copy initial conditions into the output space
  if (return_initial) {
    memcpy(write_y, y, obj->n * sizeof(double));
    store_next_output = true;
    write_y += obj->n;
  }

  if (!store_y_in_history) {
    y_next = write_y;
    out_next = write_out;
  }

  // This is used for scratch space in the case where we have output
  // (see the end of loop body).
  double *ytmp = has_output ? (double*) R_alloc(obj->n, sizeof(double)) : NULL;

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
    //
    // There is a *huge* bit of fencepost/off by one drama with
    // output; does output computed against variables y(i) go with
    // y(i) or with y(i + 1).  Hopefully this does not need to be
    // configurable because it's a real headache.
    obj->target(obj->n, obj->step, y, y_next, obj->n_out, out_next, obj->data);
    obj->step++;
    y = y_next;

    if (has_output && store_next_output) {
      if (store_y_in_history) {
        memcpy(write_out, out_next, obj->n_out * sizeof(double));
        out_next = y_next + obj->n;
      } else {
        out_next += obj->n_out;
      }
      write_out += obj->n_out;
      store_next_output = false;
    }

    if (obj->step == obj->steps[obj->steps_idx]) {
      if (store_y_in_history) {
        memcpy(write_y, y_next, obj->n * sizeof(double));
        y_next = ((double*) obj->history->head) + obj->history_idx_y;
      } else {
        y_next += obj->n;
      }
      write_y += obj->n;
      store_next_output = true;
      obj->steps_idx++;
    }

    if (store_y_in_history) {
      difeq_store_time(obj);
      double* h = (double*) ring_buffer_head_advance(obj->history);
      y_next = h + obj->history_idx_y;
      if (buffer_data != obj->history->data) {
        buffer_data = obj->history->data;
        // I think that this is always allowed because advance will
        // always move the head forward on grow, which is the only
        // case where this is triggered.
        y = y_next - obj->history_len;
        // This bit is less controversial
        out_next = y_next + obj->n_out;
      }
    }

    if (obj->step == obj->step1) {
      // The final 'y' state is used on restart
      memcpy(obj->y1, y, obj->n * sizeof(double));
      break;
    }
  }

  if (has_output && store_next_output) {
    obj->target(obj->n, obj->step, y, ytmp, obj->n_out, write_out,
                obj->data);
    if (store_y_in_history) { // NOTE: opposite direction to usual.
      memcpy(out_next, write_out, obj->n_out * sizeof(double));
    }
  }

  difeq_global_obj = NULL;
}

// These bits are all nice and don't use any globals
const double* difeq_find_step(difeq_data *obj, int step) {
  int offset = obj->step - step;
  const void *h = NULL;
  if (obj->history != NULL && offset >= 0) {
    h = ring_buffer_head_offset(obj->history, (size_t) offset);
  }
  if (h == NULL) {
    Rf_error("difeq failure: did not find step in history (at step %d)",
             (int)obj->step);
  }
  return ((double*) h) + obj->history_idx_y;
}

double yprev_1(int step, size_t i) {
  if (step <= (int) difeq_global_obj->step0) {
    return difeq_global_obj->y0[i];
  } else {
    const double * h = difeq_find_step(difeq_global_obj, step);
    // NOTE: not checking for non-null h because that's done above in
    // find_step
    return h[i];
  }
}

void yprev_all(int step, double *y) {
  if (step <= (int) difeq_global_obj->step0) {
    memcpy(y, difeq_global_obj->y0, difeq_global_obj->n * sizeof(double));
  } else {
    const double * h = difeq_find_step(difeq_global_obj, step);
    if (h != NULL) {
      memcpy(y, h, difeq_global_obj->n * sizeof(double));
    }
  }
}

void yprev_vec(int step, const size_t *idx, size_t nidx, double *y) {
  if (step <= (int) difeq_global_obj->step0) {
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

void yprev_vec_int(int step, const int *idx, size_t nidx, double *y) {
  if (step <= (int) difeq_global_obj->step0) {
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
  h[obj->history_idx_step] = obj->step;
}

void fill_na(double *x, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    x[i] = NA_REAL;
  }
}
