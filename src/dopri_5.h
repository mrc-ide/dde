#include "dopri.h"
void dopri5_step(dopri_data *obj, double h);
double dopri5_error(dopri_data *obj);
void dopri5_save_history(dopri_data *obj, double h);
double dopri5_interpolate(size_t n, double theta, double theta1,
                          const double *history);
bool dopri5_test_stiff(dopri_data *obj, double h);
