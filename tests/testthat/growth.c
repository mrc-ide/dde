#include <stddef.h>
#include <math.h>

void growth(size_t n, double t, double *y, double *dydt, void *data) {
  double *pars = (double*) data;
  double *r = pars;
  for (size_t i = 0; i < n; ++i) {
    dydt[i] = -y[i] * r[i];
  }
}

void identity(size_t n, double t, double *y, void *data) {
}

void growth_double(size_t n, double t, double *y, void *data) {
  for (size_t i = 0; i < n; ++i) {
    y[i] *= 2.0;
  }
}
