#include <dde.h>
#include <R_ext/Rdynload.h>

typedef double t_ylag_1(double, size_t);
typedef double t_ylag_all(double, double *);
typedef double t_ylag_vec(double, size_t *, size_t, double *);

double ylag_1(double t, size_t i) {
  static t_ylag_1 *fun;
  if (fun == NULL) {
    fun = (t_ylag_1*) R_GetCCallable("dde", "ylag_1");
  }
  fun(t, i);
}

double ylag_all(double t, double *y) {
  static t_ylag_all *fun;
  if (fun == NULL) {
    fun = (t_ylag_all*) R_GetCCallable("dde", "ylag_all");
  }
  fun(t, y);
}

double ylag_vec(double t, size_t *idx, size_t nidx, double *y) {
  static t_ylag_vec *fun;
  if (fun == NULL) {
    fun = (t_ylag_vec*) R_GetCCallable("dde", "ylag_vec");
  }
  fun(t, idx, n_idx, y);
}
