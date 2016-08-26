void dopri5_step(dopri5_data *obj, double h);
double dopri5_error(dopri5_data *obj);
void dopri5_save_history(dopri5_data *obj, double h);
double dopri5_interpolate(size_t n, double theta, double theta1,
                          const double *history);
