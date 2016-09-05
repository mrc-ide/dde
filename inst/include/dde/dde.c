#include <dde/dde.h>
#include <R_ext/Rdynload.h>

typedef double t_ylag_1(double, size_t);
typedef double t_ylag_all(double, double *);
typedef double t_ylag_vec(double, const size_t *, size_t, double *);
typedef double t_ylag_vec_int(double, const int *, size_t, double *);

double ylag_1(double t, size_t i) {
  static t_ylag_1 *fun;
  if (fun == NULL) {
    fun = (t_ylag_1*) R_GetCCallable("dde", "ylag_1");
  }
  return fun(t, i);
}

void ylag_all(double t, double *y) {
  static t_ylag_all *fun;
  if (fun == NULL) {
    fun = (t_ylag_all*) R_GetCCallable("dde", "ylag_all");
  }
  fun(t, y);
}

void ylag_vec(double t, const size_t *idx, size_t nidx, double *y) {
  static t_ylag_vec *fun;
  if (fun == NULL) {
    fun = (t_ylag_vec*) R_GetCCallable("dde", "ylag_vec");
  }
  fun(t, idx, nidx, y);
}

void ylag_vec_int(double t, const int *idx, size_t nidx, double *y) {
  static t_ylag_vec_int *fun;
  if (fun == NULL) {
    fun = (t_ylag_vec_int*) R_GetCCallable("dde", "ylag_vec_int");
  }
  fun(t, idx, nidx, y);
}
