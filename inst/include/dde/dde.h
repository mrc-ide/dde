#include <stddef.h>

// TODO: Consider offering a ylag() and a yprev() function (with
// trivial interfaces to t and tau respectively).
double ylag_1(double t, size_t i);
void ylag_all(double t, double *y);
void ylag_vec(double t, size_t *idx, size_t nidx, double *y);
void ylag_vec_int(double t, int *idx, size_t nidx, double *y);
