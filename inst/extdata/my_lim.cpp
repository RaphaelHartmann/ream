
#include <math.h>

/* method for the drift rate */
extern "C" double drift(const double *phi, double x, double t)  {
  double mu = phi[2];
  double l = phi[3];
  double v = mu - l*x;
  return v;
}

/* method for the diffusion rate */
extern "C" double diffusion(const double *phi, double x, double t)  {
  return phi[4];
}

/* method for the upper threshold */
extern "C" double upper_threshold(const double *phi, double t)  {
  return phi[5];
}

/* method for the lower threshold */
extern "C" double lower_threshold(const double *phi, double t)  {
  return -phi[5];
}

/* method for the contamination strength */
extern "C" double contamination_strength(const double *phi)  {
  return phi[6];
}

/* method for the contamination probability distribution */
extern "C" double contamination_probability(const double *phi, double t)  {
  double gl = phi[7];
  double gu = phi[8];
  double pg = 0.0;
  if ((t >= gl) && (t <= gu)) {
    pg = 1.0/(gu - gl);
  }
  return pg;
}
