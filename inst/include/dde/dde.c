#include <dde/dde.h>
#include <R_ext/Rdynload.h>

// dopri
typedef double t_ylag_1(double, size_t);
typedef void t_ylag_all(double, double *);
typedef void t_ylag_vec(double, const size_t *, size_t, double *);
typedef void t_ylag_vec_int(double, const int *, size_t, double *);

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

// difeq
typedef double t_yprev_1(int, size_t);
typedef void t_yprev_all(int, double *);
typedef void t_yprev_vec(int, const size_t *, size_t, double *);
typedef void t_yprev_vec_int(int, const int *, size_t, double *);

double yprev_1(int step, size_t i) {
  static t_yprev_1 *fun;
  if (fun == NULL) {
    fun = (t_yprev_1*) R_GetCCallable("dde", "yprev_1");
  }
  return fun(step, i);
}

void yprev_all(int step, double *y) {
  static t_yprev_all *fun;
  if (fun == NULL) {
    fun = (t_yprev_all*) R_GetCCallable("dde", "yprev_all");
  }
  fun(step, y);
}

void yprev_vec(int step, const size_t *idx, size_t nidx, double *y) {
  static t_yprev_vec *fun;
  if (fun == NULL) {
    fun = (t_yprev_vec*) R_GetCCallable("dde", "yprev_vec");
  }
  fun(step, idx, nidx, y);
}

void yprev_vec_int(int step, const int *idx, size_t nidx, double *y) {
  static t_yprev_vec_int *fun;
  if (fun == NULL) {
    fun = (t_yprev_vec_int*) R_GetCCallable("dde", "yprev_vec_int");
  }
  fun(step, idx, nidx, y);
}
