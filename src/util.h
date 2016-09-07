#include <R.h>
#include <Rinternals.h>

int scalar_int(SEXP x);
double scalar_double(SEXP x);
int check_index_bounds(int x, size_t len);
size_t r_index(SEXP x, size_t len);
size_t * r_indices(SEXP x, size_t len);
