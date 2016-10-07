#include <dde/dde.h>

void growth(size_t n, size_t i, double *y, double *y_next,
            size_t n_out, double *output, void *data) {
  int idx[5] = {0, 1, 2, 3, 4};
  double yprev[5];
  yprev_vec_int(i - 1, idx, 5, yprev);
  for (size_t i = 0; i < n; ++i) {
    y_next[i] = y[i] + yprev[i];
  }
}

#include <dde/dde.c>
