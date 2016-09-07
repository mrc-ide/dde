#include <stddef.h>

// dopri:
double ylag_1(double t, size_t i);
void ylag_all(double t, double *y);
void ylag_vec(double t, const size_t *idx, size_t nidx, double *y);
void ylag_vec_int(double t, const int *idx, size_t nidx, double *y);

// difeq:
double yprev_1(int step, size_t i);
void yprev_all(int step, double *y);
void yprev_vec(int step, const size_t *idx, size_t nidx, double *y);
void yprev_vec_int(int step, const int *idx, size_t nidx, double *y);
