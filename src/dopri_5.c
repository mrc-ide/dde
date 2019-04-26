#include "dopri_5.h"

// dopri5 constants
#define C2 0.2
#define C3 0.3
#define C4 0.8
#define C5 8.0 / 9.0
#define A21 0.2
#define A31 3.0 / 40.0
#define A32 9.0 / 40.0
#define A41 44.0 / 45.0
#define A42 -56.0 / 15.0
#define A43 32.0 / 9.0
#define A51 19372.0 / 6561.0
#define A52 -25360.0 / 2187.0
#define A53 64448.0 / 6561.0
#define A54 -212.0 / 729.0
#define A61 9017.0 / 3168.0
#define A62 -355.0 / 33.0
#define A63 46732.0 / 5247.0
#define A64 49.0 / 176.0
#define A65 -5103.0 / 18656.0
#define A71 35.0 / 384.0
#define A73 500.0 / 1113.0
#define A74 125.0 / 192.0
#define A75 -2187.0 / 6784.0
#define A76 11.0 / 84.0
#define E1 71.0 / 57600.0
#define E3 -71.0 / 16695.0
#define E4 71.0 / 1920.0
#define E5 -17253.0 / 339200.0
#define E6 22.0 / 525.0
#define E7 -1.0 / 40.0
// ---- DENSE OUTPUT OF SHAMPINE (1986)
#define D1 -12715105075.0 / 11282082432.0
#define D3 87487479700.0 / 32700410799.0
#define D4 -10690763975.0 / 1880347072.0
#define D5 701980252875.0 / 199316789632.0
#define D6 -1453857185.0 / 822651844.0
#define D7 69997945.0 / 29380423.0

void dopri5_step(dopri_data *obj, double h) {
  const double t = obj->t;
  const size_t n = obj->n;
  double
    *k1 = obj->k[0],
    *k2 = obj->k[1],
    *k3 = obj->k[2],
    *k4 = obj->k[3],
    *k5 = obj->k[4],
    *k6 = obj->k[5];
  double *y = obj->y, *y1 = obj->y1, *ysti = obj->k[6];

  for (size_t i = 0; i < n; ++i) { // 22
    y1[i] = y[i] + h * A21 * k1[i];
  }
  dopri_eval(obj, t + C2 * h, y1, k2);
  for (size_t i = 0; i < n; ++i) { // 23
    y1[i] = y[i] + h * (A31 * k1[i] + A32 * k2[i]);
  }
  dopri_eval(obj, t + C3 * h, y1, k3);
  for (size_t i = 0; i < n; ++i) { // 24
    y1[i] = y[i] + h * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]);
  }
  dopri_eval(obj, t + C4 * h, y1, k4);
  for (size_t i = 0; i < n; ++i) { // 25
    y1[i] = y[i] + h * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i]);
  }
  dopri_eval(obj, t + C5 * h, y1, k5);
  for (size_t i = 0; i < n; ++i) { // 26
    ysti[i] = y[i] + h * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] +
                          A64 * k4[i] + A65 * k5[i]);
  }
  double t_next = t + h;
  dopri_eval(obj, t_next, ysti, k6);
  for (size_t i = 0; i < n; ++i) { // 27
    y1[i] = y[i] + h * (A71 * k1[i] + A73 * k3[i] + A74 * k4[i] +
                        A75 * k5[i] + A76 * k6[i]);
  }
  dopri_eval(obj, t_next, y1, k2);

  // TODO: Doing this unconditionally at the moment, but this should
  // be tuned, and possibly thinned (e.g., with the index thing).
  double *history = (double*) obj->history->head;
  for (size_t i = 0, j = 4 * n; i < n; ++i, ++j) {
    history[j] =
      h * (D1 * k1[i] + D3 * k3[i] + D4 * k4[i] +
           D5 * k5[i] + D6 * k6[i] + D7 * k2[i]);
  }

  for (size_t i = 0; i < n; ++i) {
    k4[i] = h * (E1 * k1[i] + E3 * k3[i] + E4 * k4[i] +
                 E5 * k5[i] + E6 * k6[i] + E7 * k2[i]);
  }
}

double dopri5_error(dopri_data *obj) {
  double err = 0.0;
  for (size_t i = 0; i < obj->n; ++i) {
    double sk = obj->atol + obj->rtol * fmax(fabs(obj->y[i]), fabs(obj->y1[i]));
    err += square(obj->k[3][i] / sk);
  }
  return sqrt(err / obj->n);
}

void dopri5_save_history(dopri_data *obj, double h) {
  double *history = (double*) obj->history->head;
  for (size_t i = 0; i < obj->n; ++i) {
    double ydiff = obj->y1[i] - obj->y[i];
    double bspl = h * obj->k[0][i] - ydiff;
    history[             i] = obj->y[i];
    history[    obj->n + i] = ydiff;
    history[2 * obj->n + i] = bspl;
    history[3 * obj->n + i] = -h * obj->k[1][i] + ydiff - bspl;
  }
  history[obj->history_idx_time    ] = obj->t;
  history[obj->history_idx_time + 1] = h;
}

double dopri5_interpolate(size_t n, double theta, double theta1,
                          const double *history) {
  return history[0] + theta *
    (history[n] + theta1 *
     (history[2 * n] + theta *
      (history[3 * n] + theta1 *
       history[4 * n])));
}

bool dopri5_test_stiff(dopri_data *obj, double h) {
  double stnum = 0, stden = 0;
  double
    *k2 = obj->k[1],
    *k6 = obj->k[5];
  double *y1 = obj->y1, *ysti = obj->k[6];
  for (size_t i = 0; i < obj->n; ++i) {
    stnum += square(k2[i] - k6[i]);
    stden += square(y1[i] - ysti[i]);
  }
  return stden > 0 && fabs(h) * sqrt(stnum / stden) > 3.25;
}
