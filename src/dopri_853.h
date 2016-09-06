#include "dopri.h"
void dopri853_step(dopri_data *obj, double h);
double dopri853_error(dopri_data *obj);
void dopri853_save_history(dopri_data *obj, double h);
double dopri853_interpolate(size_t n, double theta, double theta1,
                            const double *history);
bool dopri853_test_stiff(dopri_data *obj, double h);
