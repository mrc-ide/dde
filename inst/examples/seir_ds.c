#include <R.h>
#include <R_ext/Rdynload.h>

void lagvalue(double tau, int *nr, int N, double *ytau);

// The parameters are going to be arranged:
//
//   t0
//   S0, I0
//   (b, N, beta, sigma, delta, t_latent)
//
// See below for why t0, S0 and I0 are stored
static double parms[3];

// The standard deSolve initialisation function
void seir_initmod(void (* odeparms)(int *, double *)) {
  int N = 3;
  odeparms(&N, parms);
}

// The RHS
void seir_deSolve(int *n, double *t, double *y, double *dydt,
                  double *yout, int *ip) {
  // again, hard-coded parameters for now; will change this shortly
  // once I get the same working with the dde impementation.
  double b = 0.1, N = 1e7, beta = 10.0, sigma = 1.0 / 3.0,
    delta = 1.0 / 21.0, t_latent = 14.0;
  double Births = N * b, surv = exp(-b * t_latent);

  // Because of the way that deSolve implements delays we need to
  // store the initial time and values in the parameters vector; if
  // the requested time is earlier than the time we started at then
  // the initial values need to be used, which we also store in the
  // parameters.
  double t0 = parms[0];
  const double tau = *t - t_latent;
  static int idx[2] = {0, 2};
  double S_lag, I_lag;
  if (tau <= t0) {
    S_lag = parms[1];
    I_lag = parms[2];
  } else {
    double ylag[2];
    lagvalue(tau, idx, 2, ylag);
    S_lag = ylag[0];
    I_lag = ylag[1];
  }

  const double S = y[0], E = y[1], I = y[2], R = y[3];
  const double new_inf = beta * S * I / N;
  const double lag_inf = beta * S_lag * I_lag * surv / N;

  dydt[0] = Births  - b * S - new_inf + delta * R;
  dydt[1] = new_inf - lag_inf - b * E;
  dydt[2] = lag_inf - (b + sigma) * I;
  dydt[3] = sigma * I - b * R - delta * R;
}

// This is the interface to deSolve's lag functions.  Note that unlike
// dde you are responsible for checking for underflows and providing
// values for underflowed times.
void lagvalue(double tau, int *nr, int N, double *ytau) {
  typedef void lagvalue_t(double, int *, int, double *);
  static lagvalue_t *fun = NULL;
  if (fun == NULL) {
    fun = (lagvalue_t*) R_GetCCallable("deSolve", "lagvalue");
  }
  fun(tau, nr, N, ytau);
}
