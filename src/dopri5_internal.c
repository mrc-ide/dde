#include "dopri5.h"
#include "util.h"
#include <R.h>

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

void dopri5_step(dopri5_data *obj, double h) {
  const double t = obj->t;
  const size_t n = obj->n;
  double
    *k1 = obj->k[0],
    *k2 = obj->k[1],
    *k3 = obj->k[2],
    *k4 = obj->k[3],
    *k5 = obj->k[4],
    *k6 = obj->k[5];

  for (size_t i = 0; i < n; ++i) { // 22
    obj->y1[i] = obj->y[i] + h * A21 * k1[i];
  }
  obj->target(n, t + C2 * h, obj->y1, k2, obj->data);
  for (size_t i = 0; i < n; ++i) { // 23
    obj->y1[i] = obj->y[i] + h * (A31 * k1[i] + A32 * k2[i]);
  }
  obj->target(n, t + C3 * h, obj->y1, k3, obj->data);
  for (size_t i = 0; i < n; ++i) { // 24
    obj->y1[i] = obj->y[i] + h * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]);
  }
  obj->target(n, t + C4 * h, obj->y1, k4, obj->data);
  for (size_t i = 0; i < n; ++i) { // 25
    obj->y1[i] = obj->y[i] + h * (A51 * k1[i] + A52 * k2[i] +
                                  A53 * k3[i] + A54 * k4[i]);
  }
  obj->target(n, t + C5 * h, obj->y1, k5, obj->data);
  for (size_t i = 0; i < n; ++i) { // 26
    obj->ysti[i] = obj->y[i] + h * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] +
                                    A64 * k4[i] + A65 * k5[i]);
  }
  double t_next = t + h;
  obj->target(n, t_next, obj->ysti, k6, obj->data);
  for (size_t i = 0; i < n; ++i) { // 27
    obj->y1[i] = obj->y[i] + h * (A71 * k1[i] + A73 * k3[i] + A74 * k4[i] +
                                  A75 * k5[i] + A76 * k6[i]);
  }
  obj->target(n, t_next, obj->y1, k2, obj->data);

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
  obj->n_eval += 6;
}
