#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

int scalar_int(SEXP x);
double scalar_double(SEXP x);
int check_index_bounds(int x, size_t len);
size_t r_index(SEXP x, size_t len);
size_t * r_indices(SEXP x, size_t len);
void * data_pointer(SEXP r_data, SEXP r_data_is_real);
void * ptr_get(SEXP r_ptr);
DL_FUNC ptr_fn_get(SEXP r_ptr);
