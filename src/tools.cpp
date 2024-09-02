/* file: model_functions.cpp
 Functions for defining model.
 Author: Mathew Murrow and Raphael Hartmann
 Date: Sep 02, 2024 */

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

#include "tools.h"


double unif_L() {

  double u;
  do {
    u = unif_rand();
  } while (u < 0 || u >= 1);
  return u;

}

/* function for the drift rate in SSP */
double ncdf(double x) {
  return 0.5 * ( 1.0 + erf(x / sqrt(2.0)) );
}
