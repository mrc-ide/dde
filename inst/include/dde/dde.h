#include <stddef.h>

// dopri:
double ylag_1(double t, size_t i);
void ylag_all(double t, double *y);
void ylag_vec(double t, const size_t *idx, size_t nidx, double *y);
void ylag_vec_int(double t, const int *idx, size_t nidx, double *y);

// difeq:
double yprev_1(size_t step, size_t i);
void yprev_all(size_t step, double *y);
void yprev_vec(size_t step, const size_t *idx, size_t nidx, double *y);
void yprev_vec_int(size_t step, const int *idx, size_t nidx, double *y);
