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

SEXP run_lorenz(SEXP r_times) {
  size_t n_times = length(r_times);
  double *times = REAL(r_times);
  double pars[3] = {10.0, 28.0, 8.0 / 3.0};
  double y[3] = {10.0, 1.0, 1.0};
  SEXP ret = PROTECT(allocVector(VECSXP, 2));
  SEXP ret_y = PROTECT(allocMatrix(REALSXP, 3, n_times - 1));
  double *y_out = REAL(ret_y);
  dopri5_data* obj = dopri5_data_alloc(&lorenz, 3, (void*) pars, 100);
  obj->rtol = 1e-7;
  obj->atol = 1e-7;
  dopri5_integrate(obj, y, times, n_times, y_out);

  Rprintf("Integration complete with code: %d (err? %d)\n",
          obj->code, obj->error);

  size_t n_history = ring_buffer_used(obj->history, 0);
  SEXP history = PROTECT(allocMatrix(REALSXP, obj->history_len, n_history));
  ring_buffer_memcpy_from(REAL(history), obj->history, n_history);
  SET_VECTOR_ELT(ret, 0, ret_y);
  SET_VECTOR_ELT(ret, 1, history);

  dopri5_data_free(obj);
  UNPROTECT(3);
  return ret;
}
