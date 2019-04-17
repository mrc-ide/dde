#include <dde/dde.h>
#include "r_dopri.h"
#include "r_difeq.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rversion.h>

static const R_CallMethodDef call_methods[] = {
  {"Cdopri",          (DL_FUNC) &r_dopri,             26},
  {"Cdopri_continue", (DL_FUNC) &r_dopri_continue,    10},
  {"Cdopri_copy",     (DL_FUNC) &r_dopri_copy,         1},
  {"Cylag",           (DL_FUNC) &r_ylag,               2},

  {"Cdifeq",          (DL_FUNC) &r_difeq,             11},
  {"Cdifeq_continue", (DL_FUNC) &r_difeq_continue,     8},
  {"Cdifeq_copy",     (DL_FUNC) &r_difeq_copy,         1},
  {"Cyprev",          (DL_FUNC) &r_yprev,              2},

  {NULL,              NULL,                            0}
};

void R_init_dde(DllInfo *info) {
  // Register C routines to be called from C:
  // dopri
  R_RegisterCCallable("dde", "ylag_1",        (DL_FUNC) &ylag_1);
  R_RegisterCCallable("dde", "ylag_all",      (DL_FUNC) &ylag_all);
  R_RegisterCCallable("dde", "ylag_vec",      (DL_FUNC) &ylag_vec);
  R_RegisterCCallable("dde", "ylag_vec_int",  (DL_FUNC) &ylag_vec_int);
  // difeq
  R_RegisterCCallable("dde", "yprev_1",       (DL_FUNC) &yprev_1);
  R_RegisterCCallable("dde", "yprev_all",     (DL_FUNC) &yprev_all);
  R_RegisterCCallable("dde", "yprev_vec",     (DL_FUNC) &yprev_vec);
  R_RegisterCCallable("dde", "yprev_vec_int", (DL_FUNC) &yprev_vec_int);

  // Register C routines to be called from R:
  R_registerRoutines(info, NULL, call_methods, NULL, NULL);

#if defined(R_VERSION) && R_VERSION >= R_Version(3, 3, 0)
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
#endif
}
