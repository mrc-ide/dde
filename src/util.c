#include "util.h"

int scalar_int(SEXP x) {
  int ret;
  if (length(x) != 1) {
    Rf_error("Expected a scalar");
  }
  if (TYPEOF(x) == INTSXP) {
    ret = INTEGER(x)[0];
  } else if (TYPEOF(x) == REALSXP) {
    ret = REAL(x)[0];
  } else {
    Rf_error("Expected an integer");
  }
  return ret;
}

double scalar_double(SEXP x) {
  double ret;
  if (length(x) != 1) {
    Rf_error("Expected a scalar");
  }
  if (TYPEOF(x) == INTSXP) {
    ret = INTEGER(x)[0];
  } else if (TYPEOF(x) == REALSXP) {
    ret = REAL(x)[0];
  } else {
    Rf_error("Expected a double");
  }
  return ret;
}

int check_index_bounds(int x, size_t len) {
  if (x < 1 || x > (int)len) {
    Rf_error("Expected index on [1..%d] (%d out of bounds)", len, x);
  }
  return x - 1;
}

size_t r_index(SEXP x, size_t len) {
  int tmp;
  if (TYPEOF(x) == INTSXP) {
    tmp = (int)INTEGER(x)[0];
  } else if (TYPEOF(x) == REALSXP) {
    tmp = (int)REAL(x)[0];
  } else {
    Rf_error("Invalid type for index");
  }
  return check_index_bounds(tmp, len);
}

size_t * r_indices(SEXP x, size_t len) {
  const size_t ni = length(x);
  size_t *idx = (size_t*) R_alloc(ni, sizeof(size_t));
  for (size_t i = 0; i < ni; ++i) {
    int tmp;
    if (TYPEOF(x) == INTSXP) {
      tmp = (int)INTEGER(x)[i];
    } else if (TYPEOF(x) == REALSXP) {
      tmp = (int)(REAL(x)[i]);
    } else {
      Rf_error("Invalid type for index");
    }
    idx[i] = check_index_bounds(tmp, len);
  }
  return idx;
}
