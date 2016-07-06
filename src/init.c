#include <dde/dde.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_dde(DllInfo *info) {
  R_RegisterCCallable("dde", "ylag_1", (DL_FUNC) &ylag_1);
  R_RegisterCCallable("dde", "ylag_all", (DL_FUNC) &ylag_all);
  R_RegisterCCallable("dde", "ylag_vec", (DL_FUNC) &ylag_vec);
  R_RegisterCCallable("dde", "ylag_vec_int", (DL_FUNC) &ylag_vec_int);
}
