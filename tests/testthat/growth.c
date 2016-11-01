#include <stddef.h>
#include <math.h>

void exponential(size_t n, double t, double *y, double *dydt, void *data) {
  double *pars = (double*) data;
  double *r = pars;
  for (size_t i = 0; i < n; ++i) {
    dydt[i] = -y[i] * r[i];
  }
}

void linear(size_t n, double t, double *y, double *dydt, void *data) {
  double *pars = (double*) data;
  double *r = pars;
    for (size_t i = 0; i < n; ++i) {
    dydt[i] = r[i];
  }
}

void identity(size_t n, double t, double *y, void *data) {
}

void double_variables(size_t n, double t, double *y, void *data) {
  for (size_t i = 0; i < n; ++i) {
    y[i] *= 2.0;
  }
}

void halve_variables(size_t n, double t, double *y, void *data) {
  for (size_t i = 0; i < n; ++i) {
    y[i] /= 2.0;
  }
}

void double_parameters(size_t n, double t, double *y, void *data) {
  double *pars = (double*) data;
  double *r = pars;
  for (size_t i = 0; i < n; ++i) {
    r[i] *= 2;
  }
}

void halve_parameters(size_t n, double t, double *y, void *data) {
  double *pars = (double*) data;
  double *r = pars;
  for (size_t i = 0; i < n; ++i) {
    r[i] /= 2;
  }
}
