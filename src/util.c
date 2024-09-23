#include "util.h"
#include <Rversion.h>

int scalar_int(SEXP x) {
  int ret;
  if (Rf_length(x) != 1) {
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
  if (Rf_length(x) != 1) {
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
    Rf_error("Expected index on [1..%d] (%d out of bounds)", (int)len, x);
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
  const size_t ni = Rf_length(x);
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

void * data_pointer(SEXP r_data, SEXP r_data_is_real) {
  void *data;
  if (TYPEOF(r_data) == REALSXP && INTEGER(r_data_is_real)[0]) {
    data = (void*) REAL(r_data);
  } else if (TYPEOF(r_data) == EXTPTRSXP) {
    data = R_ExternalPtrAddr(r_data);
    if (data == NULL) {
      // NOTE: in the underlying interface this will have come in as
      // 'parms', not 'data'.
      Rf_error("Was passed null pointer for 'parms'");
    }
  } else {
    data = (void*) r_data;
  }
  return data;
}

void * ptr_get(SEXP r_ptr) {
  if (TYPEOF(r_ptr) != EXTPTRSXP) {
    Rf_error("Expected an external pointer");
  }
  void *ret = R_ExternalPtrAddr(r_ptr);
  if (ret == NULL) {
    Rf_error("pointer has been freed (perhaps serialised?)");
  }
  return ret;
}

// This gets a function pointer from a data pointer and avoids some
// compiler warnings.
DL_FUNC ptr_fn_get(SEXP r_ptr) {
#if defined(R_VERSION) && R_VERSION >= R_Version(3, 4, 0)
  return R_ExternalPtrAddrFn(r_ptr);
#else
  return (DL_FUNC) R_ExternalPtrAddr(r_ptr);
#endif
}
