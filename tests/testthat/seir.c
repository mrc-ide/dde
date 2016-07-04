#include <dde/dde.h>
#include <stddef.h>
#include <math.h>

void seir(size_t n, double t, double *y, double *dydt, void *data) {
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

// This needs to be included exactly once per shared library.
#include <dde/dde.c>
