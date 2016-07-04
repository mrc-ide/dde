void lorenz(int n, double t, double *y, double *dydt, void *data) {
  double *pars = (double*) data;
  double sigma = pars[0];
  double R = pars[1];
  double b = pars[2];
  dydt[0] = sigma * (y[1] - y[0]);
  dydt[1] = R * y[0] - y[1] - y[0] * y[2];
  dydt[2] = -b * y[2] + y[0] * y[1];
}
