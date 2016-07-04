#include "dopri5.h"
#include <R.h>
#include <Rinternals.h>

/*
  double sq(double x) {
  return x * x;
  }

void target(int n, double t, double *y, double *dydt) {
  double amu = 0.012277471;
  double amup = 1.0 - amu;

  double r1 = sq(y[0] + amu) + sq(y[1]);
  r1 = r1 * sqrt(r1);
  double r2 = sq(y[0] - amup) + sq(y[1]);
  r2 = r2 * sqrt(r2);
  dydt[0] = y[2];
  dydt[1] = y[3];
  dydt[2] = y[0] + 2 * y[3] - amup * (y[0] + amu) / r1 -
    amu * (y[0] - amup) / r2;
  dydt[3] = y[1] - 2 * y[2] - amup * y[1] / r1 - amu * y[1] / r2;
}
*/

void lorenz(int n, double t, double *y, double *dydt, void *data) {
  double *pars = (double*) data;
  double sigma = pars[0];
  double R = pars[1];
  double b = pars[2];
  dydt[0] = sigma * (y[1] - y[0]);
  dydt[1] = R * y[0] - y[1] - y[0] * y[2];
  dydt[2] = -b * y[2] + y[0] * y[1];
}

void seir(int n, double t, double *y, double *dydt, void *data) {
  // Needs to agree with the initial conditions too; this is amazingly
  // shit (See below).  This is what we get for mimicking the deSolve
  // API though; it might be better to come up with something more
  // odin like here, where we remember what the initial conditions are
  // (especially if we're saving the initial conditions in the global
  // object this is not too bad).
  size_t I0 = 1;
  // Hard code the parameters for now, rather than passing them
  // through as `data`.
  double b = 1.0 / 10.0, N = 1e7, beta = 10.0, sigma = 1.0 / 3.0,
    delta = 1.0 / 21.0, lat_hum = 14.0;
  double Births = N * b, surv = exp(-b * lat_hum);

  // This is the same shit that deSolve does where we look up the lag
  // value from a hard coded set of numbers when we're before the
  // critical time.  It's not clever, it's something I got drawn into,
  // I'm not proud of it, but it should get the job done for now.
  //
  // The issues here are things like the requiring that t0 = 0 is
  // crap, the spray of all the duplication of initial conditions and
  // parameters, ...
  double S_lag, I_lag;
  if (t <= lat_hum) {
    S_lag = N - I0;
    I_lag = I0;
  } else {
    const double tau = t - lat_hum;
    S_lag = ylag_1(tau, 0);
    I_lag = ylag_1(tau, 2);
  }

  const double S = y[0], E = y[1], I = y[2], R = y[3];
  const double new_inf = beta * S * I / N;
  const double lag_inf = beta * S_lag * I_lag * surv / N;

  dydt[0] = Births  - b * S - new_inf + delta * R;
  dydt[1] = new_inf - lag_inf - b * E;
  dydt[2] = lag_inf - (b + sigma) * I;
  dydt[3] = sigma * I - b * R - delta * R;
}

// There are more arguments from lsoda not implemented here that will
// be needed:
//
// - rootfunc (once root-finding is supported)
// - verbose (some sort of verbose output will be useful)
// - nroot (not sure)
// - tcrit (this one is going to be needed)
// - hmin, hmax, hini
// - ynames
// - maxsteps
// - nout
// - outnames
// - forcings
// - events
//
// Some of these are big issues, some are small!

SEXP r_run_dopri5(SEXP r_y, SEXP r_times, SEXP r_func, SEXP r_data,
                  SEXP r_rtol, SEXP r_atol, SEXP data_is_real,
                  SEXP r_n_history, SEXP r_keep_history) {
  size_t n = length(r_y);
  double *y = REAL(r_y);

  size_t n_times = LENGTH(r_times);
  double *times = REAL(r_times);

  deriv_func *func = (deriv_func*)R_ExternalPtrAddr(r_func);
  void *data = INTEGER(data_is_real)[0] ? (void*) REAL(r_data) : (void*) r_data;
  size_t n_history = (size_t)INTEGER(r_n_history)[0];
  bool keep_history = INTEGER(r_keep_history)[0];

  dopri5_data* obj = dopri5_data_alloc(func, n, data, n_history);
  obj->rtol = REAL(r_rtol)[0];
  obj->atol = REAL(r_atol)[0];

  SEXP ret_y = PROTECT(allocMatrix(REALSXP, n_times - 1, n));
  dopri5_integrate(obj, y, times, n_times, REAL(ret_y));

  if (obj->error) {
    Rf_error("Integration failure with code: %d", obj->code);
  }

  if (keep_history) {
    size_t nh = ring_buffer_used(obj->history, 0);
    SEXP history = PROTECT(allocMatrix(REALSXP, obj->history_len, nh));
    ring_buffer_memcpy_from(REAL(history), obj->history, nh);
    setAttrib(ret_y, install("history"), history);
    UNPROTECT(1);
  }

  UNPROTECT(1);
  return ret_y;
}

SEXP run_lorenz(SEXP r_times, SEXP r_n_history) {
  size_t n_history = INTEGER(r_n_history)[0];
  size_t n_times = length(r_times);
  double *times = REAL(r_times);
  double pars[3] = {10.0, 28.0, 8.0 / 3.0};
  double y[3] = {10.0, 1.0, 1.0};
  SEXP ret = PROTECT(allocVector(VECSXP, 2));
  SEXP ret_y = PROTECT(allocMatrix(REALSXP, 3, n_times - 1));
  double *y_out = REAL(ret_y);
  dopri5_data* obj = dopri5_data_alloc(&lorenz, 3, (void*) pars, n_history);
  obj->rtol = 1e-7;
  obj->atol = 1e-7;
  dopri5_integrate(obj, y, times, n_times, y_out);

  Rprintf("Integration complete with code: %d (err? %d)\n",
          obj->code, obj->error);

  size_t nh = ring_buffer_used(obj->history, 0);
  SEXP history = PROTECT(allocMatrix(REALSXP, obj->history_len, nh));
  ring_buffer_memcpy_from(REAL(history), obj->history, nh);
  SET_VECTOR_ELT(ret, 0, ret_y);
  SET_VECTOR_ELT(ret, 1, history);

  dopri5_data_free(obj);
  UNPROTECT(3);
  return ret;
}
