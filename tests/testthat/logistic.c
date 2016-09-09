#include <stddef.h>
#include <math.h>

void logistic(size_t n, size_t i, double *y, double *y_next,
              size_t n_out, double *output, void *data) {
  double *pars = (double*) data;
  double *r = pars;
  for (size_t i = 0; i < n; ++i) {
    y_next[i] = r[i] * y[i] * (1 - y[i]);
  }
}
