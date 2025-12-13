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

// helper functions
double callRFunction3x(SEXP fun,
                             const double* phi, int n_phi,
                             double x, double t) {
  SEXP call, ans, phiR;
  PROTECT(phiR = Rf_allocVector(REALSXP, n_phi));
  for (int i = 0; i < n_phi; ++i)
    REAL(phiR)[i] = phi[i];
  PROTECT(call = Rf_lang4(fun, phiR,
                          Rf_ScalarReal(x),
                          Rf_ScalarReal(t)));
  PROTECT(ans = Rf_eval(call, R_GlobalEnv));   // or another environment
  double val = Rf_asReal(ans);
  UNPROTECT(3);
  return val;
}

double callRFunction2x(SEXP fun,
                             const double* phi, int n_phi, double t) {
  SEXP call, ans, phiR;
  PROTECT(phiR = Rf_allocVector(REALSXP, n_phi));
  for (int i = 0; i < n_phi; ++i)
    REAL(phiR)[i] = phi[i];
  PROTECT(call = Rf_lang3(fun, phiR, Rf_ScalarReal(t)));
  PROTECT(ans = Rf_eval(call, R_GlobalEnv));   // or another environment
  double val = Rf_asReal(ans);
  UNPROTECT(3);
  return val;
}

double callRFunction1x(SEXP fun, const double* phi, int n_phi) {
  SEXP call, ans, phiR;
  PROTECT(phiR = Rf_allocVector(REALSXP, n_phi));
  for (int i = 0; i < n_phi; ++i)
    REAL(phiR)[i] = phi[i];
  PROTECT(call = Rf_lang2(fun, phiR));
  PROTECT(ans = Rf_eval(call, R_GlobalEnv));   // or another environment
  double val = Rf_asReal(ans);
  UNPROTECT(3);
  return val;
}
