#include <R.h>
#include <Rinternals.h>

void logistic_finalize(SEXP extPtr) {
  double *p;
  if (TYPEOF(extPtr) == EXTPTRSXP) {
    p = (double*) R_ExternalPtrAddr(extPtr);
    if (p) {
      Free(p);
    }
  }
}

SEXP logistic_init(SEXP pars) {
  size_t np = length(pars);
  double *rpars = REAL(pars);
  double * p = (double*) Calloc(length(pars), double);
  memcpy(p, rpars, np * sizeof(double));

  SEXP extPtr = PROTECT(R_MakeExternalPtr(p, R_NilValue, R_NilValue));
  R_RegisterCFinalizer(extPtr, logistic_finalize);
  UNPROTECT(1);
  return extPtr;
}

void logistic(size_t n, size_t i, double *y, double *y_next,
              size_t n_out, double *output, void *data) {
  double *pars = (double*) data;
  double *r = pars;
  for (size_t i = 0; i < n; ++i) {
    y_next[i] = r[i] * y[i] * (1 - y[i]);
  }
}
