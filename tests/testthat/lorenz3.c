#include <R.h>
#include <Rinternals.h>

void lorenz_finalize(SEXP extPtr) {
  double *p;
  if (TYPEOF(extPtr) == EXTPTRSXP) {
    p = (double*) R_ExternalPtrAddr(extPtr);
    if (p) {
      Free(p);
    }
  }
}

SEXP lorenz_init(SEXP pars) {
  size_t np = length(pars);
  double *rpars = REAL(pars);
  double * p = (double*) Calloc(length(pars), double);
  memcpy(p, rpars, np * sizeof(double));

  SEXP extPtr = PROTECT(R_MakeExternalPtr(p, R_NilValue, R_NilValue));
  R_RegisterCFinalizer(extPtr, lorenz_finalize);
  UNPROTECT(1);
  return extPtr;
}

void lorenz(size_t n, double t, double *y, double *dydt, void *data) {
  double *pars = (double*) data;
  double sigma = pars[0];
  double R = pars[1];
  double b = pars[2];
  dydt[0] = sigma * (y[1] - y[0]);
  dydt[1] = R * y[0] - y[1] - y[0] * y[2];
  dydt[2] = -b * y[2] + y[0] * y[1];
}
