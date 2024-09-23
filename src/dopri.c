#include "dopri.h"
#include "dopri_5.h"
#include "dopri_853.h"
#include <R.h>
#include <stdint.h>
#include <float.h>

dopri_data* dopri_data_alloc(deriv_func* target, size_t n,
                             output_func* output, size_t n_out,
                             void *data,
                             dopri_method method,
                             size_t n_history, bool grow_history,
                             dopri_verbose verbose, SEXP  callback) {
  dopri_data *ret = (dopri_data*) R_Calloc(1, dopri_data);
  overflow_action on_overflow =
    grow_history ? OVERFLOW_GROW : OVERFLOW_OVERWRITE;

  ret->target = target;
  ret->output = output;
  ret->data = data;

  ret->method = method;
  ret->order = ret->method == DOPRI_5 ? 5 : 8;

  ret->n = n;
  ret->n_out = n_out;

  ret->n_times = 0;
  ret->times = NULL;
  ret->tcrit = NULL;

  ret->verbose = verbose;
  ret->callback = callback;

  // tcrit variables are set in reset

  // State vectors
  ret->y0 = R_Calloc(n, double); // initial
  ret->y  = R_Calloc(n, double); // current
  ret->y1 = R_Calloc(n, double); // next

  // NOTE: There's no real reason to believe that the storage
  // requirements (nk) will always grow linearly like this, but I
  // don't really anticipate adding any other schemes soon anyway, so
  // the fact that this works well for the two we have is enough.
  size_t nk = ret->order + 2;
  ret->k = R_Calloc(nk, double*);
  for (size_t i = 0; i < nk; ++i) {
    ret->k[i] = R_Calloc(n, double);
  }

  ret->history_len = 2 + ret->order * n;
  ret->history = ring_buffer_create(n_history,
                                    ret->history_len * sizeof(double),
                                    on_overflow);
  ret->history_idx_time = ret->order * n;

  // NOTE: The numbers below are defaults only, but to alter them,
  // directly modify the object after creation.  I may set up some
  // sort of helper functions for this later, but for now just don't
  // set anything stupid.
  //
  // TODO: Support vectorised tolerances?
  ret->atol = 1e-6;
  ret->rtol = 1e-6;

  switch (ret->method) {
  case DOPRI_5:
    ret->step_factor_min = 0.2;  // from dopri5.f:276, retard.f:328
    ret->step_factor_max = 10.0; // from dopri5.f:281, retard.f:333
    ret->step_beta = 0.04;       // from dopri5.f:287, retard.f:339
    break;
  case DOPRI_853:
    ret->step_factor_min = 0.333; // from dopri853.f:285
    ret->step_factor_max = 6.0;   // from dopri853.f:290
    ret->step_beta = 0.0;         // from dopri853.f:296
    break;
  }
  ret->step_size_min = DBL_EPSILON;
  ret->step_size_max = DBL_MAX;
  ret->step_size_initial = 0.0;
  ret->step_size_min_allow = false;
  ret->step_max_n = 100000;    // from dopri5.f:212
  ret->step_factor_safe = 0.9; // from dopri5.f:265

  // How often to check for stiffness?
  ret->stiff_check = 0;

  return ret;
}

dopri_data* dopri_data_copy(const dopri_data* obj) {
  size_t n_history = ring_buffer_size(obj->history, false);
  bool grow_history = obj->history->on_overflow == OVERFLOW_GROW;
  dopri_data* ret = dopri_data_alloc(obj->target, obj->n,
                                     obj->output, obj->n_out,
                                     obj->data, obj->method,
                                     n_history, grow_history,
                                     obj->verbose, obj->callback);
  // Then update a few things
  ret->t0 = obj->t0;
  ret->t  = obj->t;
  ring_buffer_mirror(obj->history, ret->history);
  ret->history_idx_time = obj->history_idx_time;
  ret->sign = obj->sign;
  ret->atol = obj->atol;
  ret->rtol = obj->rtol;
  ret->step_factor_safe = obj->step_factor_safe;
  ret->step_factor_min = obj->step_factor_min;
  ret->step_factor_max = obj->step_factor_max;
  ret->step_size_min = obj->step_size_min;
  ret->step_size_max = obj->step_size_max;
  ret->step_size_initial = obj->step_size_initial;
  ret->step_max_n = obj->step_max_n;
  ret->step_size_min_allow = obj->step_size_min_allow;
  ret->step_beta = obj->step_beta;
  ret->step_constant = obj->step_constant;

  // NOTE: I'm not copying the times over because the first thing that
  // happens is always a reset and that deals with this.
  ret->times = NULL;
  ret->tcrit = NULL;

  // Copy all of the internal state:
  memcpy(ret->y, obj->y, obj->n * sizeof(double));
  memcpy(ret->y0, obj->y0, obj->n * sizeof(double));
  memcpy(ret->y1, obj->y1, obj->n * sizeof(double));
  size_t nk = ret->order + 2;
  for (size_t i = 0; i < nk; ++i) {
    memcpy(ret->k[i], obj->k[i], obj->n * sizeof(double));
  }

  return ret;
}

// We'll need a different reset when we're providing history, because
// then we won't end up resetting t0/y0 the same way.
void dopri_data_reset(dopri_data *obj, const double *y,
                      const double *times, size_t n_times,
                      const double *tcrit, size_t n_tcrit,
                      const bool *is_event, event_func **events) {
  obj->error = false;
  obj->code = NOT_SET;

  switch (obj->method) {
  case DOPRI_5:
    obj->step_constant = 0.2 - obj->step_beta * 0.75;
    break;
  case DOPRI_853:
    // NOTE: probably worth chasing down the source of these
    // calculations.  It does look like 1/order - beta * C, but not
    // sure where C comes from.
    obj->step_constant = 0.125 - obj->step_beta * 0.2;
    break;
  }

  memcpy(obj->y0, y, obj->n * sizeof(double));
  if (obj->y != y) { // this is true on some restarts
    memcpy(obj->y, y, obj->n * sizeof(double));
  }

  obj->n_times = n_times;
  obj->times = times;
  obj->times_idx = 1; // skipping the first time!

  // TODO: I don't check that there is at least one time anywhere in
  // *this* routine, but it is checked in r_dopri which calls this.
  if (times[n_times - 1] == times[0]) {
    obj->error = true;
    obj->code = ERR_ZERO_TIME_DIFFERENCE;
    return;
  }
  obj->sign = copysign(1.0, times[n_times - 1] - times[0]);
  for (size_t i = 0; i < n_times - 1; ++i) {
    double t0 = times[i], t1 = times[i + 1];
    bool err = (obj->sign > 0 && t1 < t0) || (obj->sign < 0 && t1 > t0);
    // perhaps (obj->sign) * (t1 < t0) < 0
    if (err) {
      obj->error = true;
      obj->code = ERR_INCONSISTENT_TIME;
      return;
    }
  }

  if (ring_buffer_is_empty(obj->history)) {
    obj->t0 = obj->sign * times[0];
  } else {
    // This is the restart condition
    double * h = (double*) ring_buffer_tail(obj->history);
    obj->t0 = obj->sign * h[obj->history_idx_time];
  }
  obj->t = times[0];

  obj->n_tcrit = n_tcrit;
  obj->tcrit = tcrit;
  obj->tcrit_idx = 0;
  if (n_tcrit > 0) {
    double t0 = obj->sign * times[0]; // because of the restart condition above.
    while (obj->tcrit_idx < n_tcrit &&
           obj->sign * tcrit[obj->tcrit_idx] <= t0) {
      obj->tcrit_idx++;
    }
  }

  obj->is_event = is_event;
  obj->events = events;

  obj->n_eval = 0;
  obj->n_step = 0;
  obj->n_accept = 0;
  obj->n_reject = 0;

  if (obj->stiff_check == 0) {
    obj->stiff_check = SIZE_MAX;
  }

  obj->stiff_n_stiff = 0;
  obj->stiff_n_nonstiff = 0;
}

void dopri_data_free(dopri_data *obj) {
  R_Free(obj->y0);
  R_Free(obj->y);
  R_Free(obj->y1);

  size_t nk = obj->order + 2;
  for (size_t i = 0; i < nk; ++i) {
    R_Free(obj->k[i]);
  }
  R_Free(obj->k);

  ring_buffer_destroy(obj->history);

  R_Free(obj);
}

// This is super ugly, but needs to be done so that the lag functions
// can access the previous history easily.  I don't see an obvious way
// around this unfortunately, given that the lag functions need to be
// callable in user code (so without forcing some weird and blind
// passing of a void object around, which will break compatibility
// with deSolve even more and make the interface for the dde and
// non-dde equations quite different) this seems like a reasonable way
// of achiving this.  Might change later though.
static dopri_data *dopri_global_obj;

// Used to query the problem size safely from the r_dopri.c file
size_t get_current_problem_size_dde(void) {
  return dopri_global_obj == NULL ? 0 : dopri_global_obj->n;
}

// Wrappers around the two methods:
void dopri_step(dopri_data *obj, double h) {
  switch (obj->method) {
  case DOPRI_5:
    dopri5_step(obj, h);
    break;
  case DOPRI_853:
    dopri853_step(obj, h);
    break;
  }
  dopri_print_step(obj, h);
}

double dopri_error(dopri_data *obj) {
  double ret = 0; // avoids a compiler warning
  switch (obj->method) {
  case DOPRI_5:
    ret = dopri5_error(obj);
    break;
  case DOPRI_853:
    ret = dopri853_error(obj);
    break;
  }
  return ret;
}

void dopri_save_history(dopri_data *obj, double h) {
  switch (obj->method) {
  case DOPRI_5:
    dopri5_save_history(obj, h);
    break;
  case DOPRI_853:
    dopri853_save_history(obj, h);
    break;
  }
}

// Integration is going to be over a set of times 't', of which there
// are 'n_t'.
void dopri_integrate(dopri_data *obj, const double *y,
                     const double *times, size_t n_times,
                     const double *tcrit, size_t n_tcrit,
                     const bool *is_event, event_func **events,
                     double *y_out, double *out,
                     bool return_initial) {
  dopri_data_reset(obj, y, times, n_times, tcrit, n_tcrit,
                   is_event, events);
  if (obj->error) {
    return;
  }

  double fac_old = 1e-4;
  bool stop = false, last = false, reject = false;

  double t_end = times[n_times - 1];
  double t_stop = t_end;
  if (obj->tcrit_idx < obj->n_tcrit &&
      obj->sign * obj->tcrit[obj->tcrit_idx] < obj->sign * t_end) {
    t_stop = obj->tcrit[obj->tcrit_idx];
  } else {
    t_stop = t_end;
  }

  // Possibly only set this if the number of history variables is
  // nonzero?  Needs to be set before any calls to target() though.
  dopri_global_obj = obj;

  // If requested, copy initial conditions into the output space
  if (return_initial) {
    memcpy(y_out, y, obj->n * sizeof(double));
    y_out += obj->n;
    if (obj->n_out > 0) {
      obj->output(obj->n, times[0], y, obj->n_out, out, obj->data);
      out += obj->n_out;
    }
  }

  dopri_eval(obj, obj->t, obj->y, obj->k[0]); // y => k1

  for (size_t i = 0; i < obj->n; ++i) {
    if (!R_FINITE(obj->k[0][i])) {
      Rf_error("non-finite derivative at initial time for element %d\n",
               (int)i + 1);
    }
  }

  // Work out the initial step size:
  double h = dopri_h_init(obj);
  double h_save = 0.0;

  // TODO: factor into its own thing
  while (true) {
    bool force_this_step = false;
    if (obj->n_step > obj->step_max_n) {
      obj->error = true;
      obj->code = ERR_TOO_MANY_STEPS;
      break;
    }
    if (fabs(h) <= obj->step_size_min) {
      if (obj->step_size_min_allow) {
        h = copysign(obj->step_size_min, obj->sign);
        force_this_step = true;
      } else {
        obj->error = true;
        obj->code = ERR_STEP_SIZE_TOO_SMALL;
        break;
      }
    } else if (fabs(h) <= fabs(obj->t) * DBL_EPSILON) {
      obj->error = true;
      obj->code = ERR_STEP_SIZE_VANISHED;
      break;
    }
    if ((obj->t + 1.01 * h - t_stop) * obj->sign > 0.0) {
      h = t_stop - obj->t;
      stop = true;
    }
    if ((obj->t + 1.01 * h - t_end) * obj->sign > 0.0) {
      h_save = h;
      h = t_end - obj->t;
      last = true;
    }
    // NOTE: in retard.f there is an else condition:
    //
    //   else if ((t + 1.8 * h - t_end) * sign > 0) {
    //     h = (t_end - t) * 0.55
    //   }
    //
    obj->n_step++;

    // TODO: In the Fortran there is an option here to check the irtrn
    // flag for a code of '2' which indicates that the variables need
    // recomputing.  This would be the case possibly where output
    // variables are being calculated so might want to put this back
    // in once the signalling is done.
    dopri_step(obj, h);
    if (obj->error) {
      break;
    }

    // Error estimation:
    double err = dopri_error(obj);
    double h_new = dopri_h_new(obj, fac_old, h, err);
    const bool accept = err <= 1;

    if (accept || force_this_step) {
      // Step is accepted :)
      fac_old = fmax(err, 1e-4);
      obj->n_accept++;
      if (obj->method == DOPRI_853) {
        double *k4 = obj->k[3], *k5 = obj->k[4];
        dopri_eval(obj, obj->t, k5, k4);
      }
      if (dopri_test_stiff(obj, h)) {
        obj->error = true;
        obj->code = ERR_STIFF;
        return;
      }

      dopri_save_history(obj, h);

      // TODO: it's quite possible that we can swap the pointers here
      // and avoid the memcpy.
      switch (obj->method) {
      case DOPRI_5:
        memcpy(obj->k[0], obj->k[1], obj->n * sizeof(double)); // k1 = k2
        memcpy(obj->y,    obj->y1,   obj->n * sizeof(double)); // y  = y1
        break;
      case DOPRI_853:
        memcpy(obj->k[0], obj->k[3], obj->n * sizeof(double)); // k1 = k4
        memcpy(obj->y,    obj->k[4], obj->n * sizeof(double)); // y  = k5
        break;
      }
      obj->t += h;

      // The next six lines contain a workaround for a problem in which
      // gcc 4.9.3 on Win32 incorrectly optimises this (for -O1, -O2, and -O3)
      // changing the behaviour and causing too few iterations to take place
      // under certain circumstances.

      // The use of the variable `last`, along with the way the while loop
      // logic is now written, together seem to persuade gcc against whatever
      // optimisation causes the problem; it is difficult to tell exactly how.

      // See https://github.com/mrc-ide/dde/issues/14 for the original issue
      // and https://github.com/mrc-ide/dde/pull/19 for the specific changes.
      while (last ||
             obj->sign * obj->times[obj->times_idx] <= obj->sign * obj->t) {
        if (obj->times_idx >= obj->n_times) {
          // Exists so that we eventually exit on the 'last' integration step.
          break;
        }
        // Here, it might be nice to allow transposed output or not;
        // that would be an argument to interpolate_all.  That's a bit
        // of a faff.
        dopri_interpolate_all((double *) obj->history->head, obj->method,
                               obj->n, obj->times[obj->times_idx], y_out);
        if (obj->n_out > 0) {
          obj->output(obj->n, obj->times[obj->times_idx], y_out,
                      obj->n_out, out, obj->data);
          out += obj->n_out;
        }

        y_out += obj->n;
        obj->times_idx++;
      }

      // Advance the ring buffer; we'll write to the next place after
      // this.
      ring_buffer_head_advance(obj->history);

      if (last) {
        obj->step_size_initial = h_save;
        obj->code = OK_COMPLETE;
        break;
      }
      // TODO: To understand this bit I think I will need to get the
      // book and actually look at the dopri integration bit.
      if (fabs(h_new) >= obj->step_size_max) {
        h_new = copysign(obj->step_size_max, obj->sign);
      }
      if (reject) {
        h_new = copysign(fmin(fabs(h_new), fabs(h)), obj->sign);
        reject = false;
      }
      if (stop) {
        while (obj->tcrit_idx < n_tcrit &&
               obj->sign * tcrit[obj->tcrit_idx] <= obj->sign * obj->t) {
          if (obj->is_event[obj->tcrit_idx]) {
            event_func *event = obj->events[obj->tcrit_idx];
            event(obj->n, obj->t, obj->y, obj->data);
          }
          obj->tcrit_idx++;
        }
        if (obj->tcrit_idx < obj->n_tcrit &&
            obj->sign * obj->tcrit[obj->tcrit_idx] < obj->sign * t_end) {
          t_stop = obj->tcrit[obj->tcrit_idx];
        } else {
          t_stop = t_end;
        }
        stop = false;
      } else {
        h = h_new;
      }
    } else {
      // Step is rejected :(
      //
      // TODO: This is very annoying because we need to re-invert the
      // minimum step factor.
      double fac11 = pow(err, obj->step_constant);
      h_new = h / fmin(1 / obj->step_factor_min, fac11 / obj->step_factor_safe);
      reject = true;
      if (obj->n_accept >= 1) {
        obj->n_reject++;
      }
      last = false;
      stop = false;
      h = h_new;
    }
  }

  // Reset the global state
  dopri_global_obj = NULL;
}

double dopri_h_new(dopri_data *obj, double fac_old, double h, double err) {
  double fac11 = pow(err, obj->step_constant);
  double step_factor_min = 1.0 / obj->step_factor_min;
  double step_factor_max = 1.0 / obj->step_factor_max;
  // Lund-stabilisation
  double fac = fac11 / pow(fac_old, obj->step_beta);
  fac = fmax(step_factor_max,
             fmin(step_factor_min, fac / obj->step_factor_safe));
  return h / fac;
}

double dopri_h_init(dopri_data *obj) {
  if (obj->step_size_initial > 0.0) {
    return obj->step_size_initial;
  }

  // NOTE: This is destructive with respect to most of the information
  // in the object; in particular k2, k3 will be modified.
  double *f0 = obj->k[0];

  // Compute a first guess for explicit Euler as
  //   h = 0.01 * norm (y0) / norm (f0)
  // the increment for explicit euler is small compared to the solution
  double norm_f = 0.0, norm_y = 0.0;
  for (size_t i = 0; i < obj->n; ++i) {
    double sk = obj->atol + obj->rtol * fabs(obj->y[i]);
    norm_f += square(f0[i] / sk);
    norm_y += square(obj->y[i]  / sk);
  }

  double h = (norm_f <= 1e-10 || norm_y <= 1e-10) ?
    1e-6 : sqrt(norm_y / norm_f) * 0.01;
  h = copysign(fmin(h, obj->step_size_max), obj->sign);

  double *f1 = obj->k[1], *y1 = obj->k[2];
  // Perform an explicit Euler step
  for (size_t i = 0; i < obj->n; ++i) {
    y1[i] = obj->y[i] + h * f0[i];
  }
  dopri_eval(obj, obj->t + h, y1, f1);

  // Estimate the second derivative of the solution:
  double der2 = 0.0;
  for (size_t i = 0; i < obj->n; ++i) {
    double sk = obj->atol + obj->rtol * fabs(obj->y[i]);
    der2 += square((f1[i] - f0[i]) / sk);
  }
  der2 = sqrt(der2) / h;

  // Step size is computed such that
  //   h^order * fmax(norm(f0), norm(der2)) = 0.01
  double der12 = fmax(fabs(der2), sqrt(norm_f));
  double h1 = (der12 <= 1e-15) ?
    fmax(1e-6, fabs(h) * 1e-3) : pow(0.01 / der12, 1.0 / obj->order);
  h = fmin(fmin(100 * fabs(h), h1), obj->step_size_max);
  return copysign(h, obj->sign);
}

// Specific...

// There are several interpolation functions here;
//
// * interpolate_1; interpolate a single variable i
// * interpolate_all: interpolate the entire vector
// * interpolate_idx: interpolate some of the vector
// * interpolate_idx_int: As for _idx but with an integer index (see below)
double dopri_interpolate_1(const double *history, dopri_method method,
                           size_t n, double t, size_t i) {
  const size_t idx_t = (method == DOPRI_5 ? 5 : 8) * n;
  const double t_old = history[idx_t], h = history[idx_t + 1];
  const double theta = (t - t_old) / h;
  const double theta1 = 1 - theta;
  double ret = 0.0;
  switch (method) {
  case DOPRI_5:
    ret = dopri5_interpolate(n, theta, theta1, history + i);
    break;
  case DOPRI_853:
    ret = dopri853_interpolate(n, theta, theta1, history + i);
    break;
  }
  return ret;
}

void dopri_interpolate_all(const double *history, dopri_method method,
                           size_t n, double t, double *y) {
  const size_t idx_t = (method == DOPRI_5 ? 5 : 8) * n;
  const double t_old = history[idx_t], h = history[idx_t + 1];
  const double theta = (t - t_old) / h;
  const double theta1 = 1 - theta;
  switch (method) {
  case DOPRI_5:
    for (size_t i = 0; i < n; ++i) {
      y[i] = dopri5_interpolate(n, theta, theta1, history + i);
    }
    break;
  case DOPRI_853:
    for (size_t i = 0; i < n; ++i) {
      y[i] = dopri853_interpolate(n, theta, theta1, history + i);
    }
    break;
  }
}

void dopri_interpolate_idx(const double *history, dopri_method method,
                           size_t n, double t, const size_t * idx, size_t nidx,
                           double *y) {
  const size_t idx_t = (method == DOPRI_5 ? 5 : 8) * n;
  const double t_old = history[idx_t], h = history[idx_t + 1];
  const double theta = (t - t_old) / h;
  const double theta1 = 1 - theta;
  switch (method) {
  case DOPRI_5:
    for (size_t i = 0; i < nidx; ++i) {
      y[i] = dopri5_interpolate(n, theta, theta1, history + idx[i]);
    }
    break;
  case DOPRI_853:
    for (size_t i = 0; i < nidx; ++i) {
      y[i] = dopri853_interpolate(n, theta, theta1, history + idx[i]);
    }
    break;
  }
}

// This exists to deal with deSolve taking integer arguments (and
// therefore messing up odin).  I could rewrite the whole thing here
// to use int, but that seems needlessly crap.  The issue here is only
// the *pointer* '*idx' and not anything else because we can safely
// cast the plain data arguments.  This affects only this function as
// it's the only one that takes size_t*
void dopri_interpolate_idx_int(const double *history, dopri_method method,
                               size_t n, double t, const int *idx, size_t nidx,
                               double *y) {
  const size_t idx_t = (method == DOPRI_5 ? 5 : 8) * n;
  const double t_old = history[idx_t], h = history[idx_t + 1];
  const double theta = (t - t_old) / h;
  const double theta1 = 1 - theta;
  switch (method) {
  case DOPRI_5:
    for (size_t i = 0; i < nidx; ++i) {
      y[i] = dopri5_interpolate(n, theta, theta1, history + idx[i]);
    }
    break;
  case DOPRI_853:
    for (size_t i = 0; i < nidx; ++i) {
      y[i] = dopri853_interpolate(n, theta, theta1, history + idx[i]);
    }
    break;
  }
}

bool dopri_test_stiff(dopri_data *obj, double h) {
  bool ret = false;
  if (obj->stiff_n_stiff > 0 || obj->n_accept % obj->stiff_check == 0) {
    bool stiff = false;
    switch (obj->method) {
    case DOPRI_5:
      stiff = dopri5_test_stiff(obj, h);
      break;
    case DOPRI_853:
      stiff = dopri853_test_stiff(obj, h);
      break;
    }

    if (stiff) {
      obj->stiff_n_nonstiff = 0;
      if (obj->stiff_n_stiff++ >= 15) {
        ret = true;
      }
    } else if (obj->stiff_n_stiff > 0) {
      if (obj->stiff_n_nonstiff++ >= 6) {
        obj->stiff_n_stiff = 0;
      }
    }
  }
  return ret;
}

// history searching:
//
// The big challenge here is that the history object will need to be
// found.  deSolve deals with this generally by stashing a copy of the
// objects as a global, This is nice because then the target function
// does not need to know about the struct; if you require that the
// struct is present then the whole thing falls apart a little bit
// because the target function is in the main struct, but needs to
// *call* the main struct to get the history.  So, this is a bit
// tricky to get right.

// The global/lock solver approach of deSolve is reasonable, and while
// inelegant should serve us reasonably well.  The other approach
// would be to somehow do a more C++ approach where things know more
// about the solution but that seems shit.

// So I'll need to get this right.  On integration we'll have a
// `global_state` variable that is (possibly static) initialised when
// the integration starts, and which holds a pointer to the
// integration data.  Then when the integration runs any call to ylag*
// will end up triggering that.  On exit of the itegrator we'll set
// that to NULL and hope that there are no longjup style exists
// anywhere so we can use `==NULL` safely.

// TODO: There's an issue here where if the chosen time is the last
// time (a tail offset of n-1) then we might actually need to
// interpolate off of the head which is not actually returned.  This
// is a bit of a shit, and might be worth treating separately in the
// search, or in the code below.

// These bits are all nice and don't use any globals
struct dopri_find_time_pred_data {
  size_t idx;
  double time;
};

bool dopri_find_time_forward(const void *x, void *data) {
  const struct dopri_find_time_pred_data *d =
    (struct dopri_find_time_pred_data *) data;
  return ((double*) x)[d->idx] <= d->time;
}
bool dopri_find_time_backward(const void *x, void *data) {
  const struct dopri_find_time_pred_data *d =
    (struct dopri_find_time_pred_data *) data;
  return ((double*) x)[d->idx] >= d->time;
}

const double* dopri_find_time(dopri_data *obj, double t) {
  const size_t idx_t = obj->history_idx_time;
  struct dopri_find_time_pred_data data = {idx_t, t};
  // The first shot at idx here is based on a linear interpolation of
  // the time; hopefully this gets is close to the correct point
  // without having to have a really long search time.
  const size_t n = ring_buffer_used(obj->history, 0);
  size_t idx0 = 0;
  if (n > 0) {
    const double
      t0 = ((double*) ring_buffer_tail(obj->history))[idx_t],
      t1 = ((double*) ring_buffer_tail_offset(obj->history, n - 1))[idx_t];
    if ((t0 - t) * (t1 - t) < 0) {
      idx0 = (t - t0) / (t1 - t0) / (n - 1);
    }
  }
  const void *h =
    ring_buffer_search_bisect(obj->history, idx0,
                              obj->sign > 0 ?
                              &dopri_find_time_forward :
                              &dopri_find_time_backward,
                              &data);
  if (h == NULL) {
    obj->error = true;
    obj->code = ERR_YLAG_FAIL;
  }
  return (double*) h;
}

// But these all use the global state object (otherwise these all pick
// up a `void *data` argument which will be cast to `dopri_data*`,
// but then the derivative function needs the same thing, which is
// going to seem weird and also means that the same function can't be
// easily used for dde and non dde use).
double ylag_1(double t, size_t i) {
  if (dopri_global_obj->sign * t <= dopri_global_obj->t0) {
    return dopri_global_obj->y0[i];
  } else {
    const double * h = dopri_find_time(dopri_global_obj, t);
    if (h == NULL) {
      return NA_REAL;
    } else {
      return dopri_interpolate_1(h, dopri_global_obj->method,
                                 dopri_global_obj->n, t, i);
    }
  }
}

void ylag_all(double t, double *y) {
  if (dopri_global_obj->sign * t <= dopri_global_obj->t0) {
    memcpy(y, dopri_global_obj->y0, dopri_global_obj->n * sizeof(double));
  } else {
    const double * h = dopri_find_time(dopri_global_obj, t);
    if (h != NULL) {
      dopri_interpolate_all(h, dopri_global_obj->method,
                            dopri_global_obj->n, t, y);
    }
  }
}

void ylag_vec(double t, const size_t *idx, size_t nidx, double *y) {
  if (dopri_global_obj->sign * t <= dopri_global_obj->t0) {
    for (size_t i = 0; i < nidx; ++i) {
      y[i] = dopri_global_obj->y0[idx[i]];
    }
  } else {
    const double * h = dopri_find_time(dopri_global_obj, t);
    if (h != NULL) {
      dopri_interpolate_idx(h, dopri_global_obj->method,
                            dopri_global_obj->n, t, idx, nidx, y);
    }
  }
}

void ylag_vec_int(double t, const int *idx, size_t nidx, double *y) {
  if (dopri_global_obj->sign * t <= dopri_global_obj->t0) {
    for (size_t i = 0; i < nidx; ++i) {
      y[i] = dopri_global_obj->y0[idx[i]];
    }
  } else {
    const double * h = dopri_find_time(dopri_global_obj, t);
    if (h != NULL) {
      dopri_interpolate_idx_int(h, dopri_global_obj->method,
                                dopri_global_obj->n, t, idx, nidx, y);
    }
  }
}


void dopri_eval(dopri_data *obj, double t, double *y, double *dydt) {
  dopri_print_eval(obj, t, y);
  obj->target(obj->n, t, y, dydt, obj->data);
  obj->n_eval++;
}


void dopri_print_step(dopri_data *obj, double h) {
  if (obj->verbose >= VERBOSE_STEP) {
    if (obj->callback == R_NilValue) {
      Rprintf("[step] t: %f, h: %e\n", obj->t, h);
    } else {
      dopri_callback(obj, obj->t, h, obj->y);
    }
  }
}


void dopri_print_eval(dopri_data *obj, double t, double *y) {
  if (obj->verbose >= VERBOSE_EVAL) {
    if (obj->callback == R_NilValue) {
      Rprintf("[eval] t: %f\n", t);
    } else {
      dopri_callback(obj, t, NA_REAL, y);
    }
  }
}


void dopri_callback(dopri_data *obj, double t, double h, double *y) {
  SEXP callback = VECTOR_ELT(obj->callback, 0);
  SEXP env = VECTOR_ELT(obj->callback, 1);

  SEXP r_t = PROTECT(Rf_ScalarReal(t));
  SEXP r_h = PROTECT(Rf_ScalarReal(h));
  SEXP r_y = PROTECT(Rf_allocVector(REALSXP, obj->n));
  memcpy(REAL(r_y), y, obj->n * sizeof(double));
  SEXP call = PROTECT(Rf_lang4(callback, r_t, r_h, r_y));
  Rf_eval(call, env);
  UNPROTECT(4);
}


// Utility
double square(double x) {
  return x * x;
}
