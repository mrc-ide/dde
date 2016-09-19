#include <dde/dde.h>
#include <stddef.h>
#include <math.h>

void seir(size_t n, double t, double *y, double *dydt, void *data) {
  // Hard code the parameters for now, rather than passing them
  // through as `data`.
  double b = 1.0 / 10.0, N = 1e7, beta = 10.0, sigma = 1.0 / 3.0,
    delta = 1.0 / 21.0, lat_hum = 14.0;
  double Births = N * b, surv = exp(-b * lat_hum);

  const double tau = t - lat_hum;

  static const int idx[2] = {0, 2};
  double ylag[2];
  ylag_vec_int(tau, idx, 2, ylag);
  double S_lag = ylag[0], I_lag = ylag[1];

  // This is an alternative mode that looks up the value of the lags
  // by index.
  //   size_t lag_idx[2] = {0, 2};
  //   double y_lag[2];
  //   ylag_vec(tau, lag_idx, 2, y_lag);
  //   double S_lag = y_lag[0], I_lag = y_lag[1];

  const double S = y[0], E = y[1], I = y[2], R = y[3];
  const double new_inf = beta * S * I / N;
  const double lag_inf = beta * S_lag * I_lag * surv / N;

  dydt[0] = Births  - b * S - new_inf + delta * R;
  dydt[1] = new_inf - lag_inf - b * E;
  dydt[2] = lag_inf - (b + sigma) * I;
  dydt[3] = sigma * I - b * R - delta * R;
}

// An output function that takes the sum over all compartments in the
// model.  Unlike deSolve, these can be added and removed at runtime,
// rather than just at compile time.
//
// The downside of this is that it's going to require an additional
// set of calls to the lag functions to get access to lag variables.
// That's going to complicate things a bit for odin I think.
void seir_output(size_t n, double t, const double *y,
                 size_t n_out, double *out, const void *data) {
  out[0] = y[0] + y[1] + y[2] + y[3];
}

// This needs to be included exactly once per shared library.
#include <dde/dde.c>
