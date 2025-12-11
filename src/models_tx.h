/* file: model_functions.cpp
 Functions for defining model.
 Author: Mathew Murrow and Raphael Hartmann
 Date: Sep 02, 2024 */

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

#ifndef MODELS_TX_H
#define MODELS_TX_H

#include "Model_TX.h"
#include <Rinternals.h>


class LIMF : public Model_TX {
protected:

  /* method for the non-decision time */
  double non_decision(const double phi[11]) const override {
    return phi[0];
  }

  /* method for the start point */
  double relative_start(const double phi[11]) const override {
    return phi[1];
  }

  /* method for the drift rate */
  double drift(const double phi[11], double x, double t) const override {
    double mu1 = phi[2];
    double mu2 = phi[3];
    double l = phi[4];
    double t0 = phi[5];
    double v = 0.0;

    if (t < t0) {
      v = mu1 - l*x;
    } else {
      v = mu2 - l*x;
    }

    return v;
  }

  /* method for the diffusion rate */
  double diffusion(const double phi[11], double x, double t) const override {
    return phi[6];
  }

  /* method for the upper threshold */
  double upper_threshold(const double phi[11], double t) const override {
    return phi[7];
  }

  /* method for the lower threshold */
  double lower_threshold(const double phi[11], double t) const override {
    return -phi[7];
  }

  /* method for the contamination strength */
  double contamination_strength(const double phi[11]) const override {
    return phi[8];
  }

  /* method for the contamination probability distribution */
  double contamination_probability(const double phi[11], double t) const override {
    double gl = phi[9];
    double gu = phi[10];
    double pg = 0.0;
    if ((t >= gl) && (t <= gu)) {
      pg = 1.0/(gu - gl);
    }
    return pg;
  }

  /* method for locally modifying the time step size */
  double modify_dt(const double phi[11], double t) const override {
    double t0 = phi[5];
    double range = 0.1;
    double dt_mod = 1.0;

    if ((t >= t0 - range) && (t <= t0 + range)) {
      dt_mod = 0.1;
    }

    return dt_mod;
  }

};



class LIM : public Model_TX {
protected:

  /* method for the non-decision time */
  double non_decision(const double phi[9]) const override {
    return phi[0];
  }

  /* method for the start point */
  double relative_start(const double phi[9]) const override {
    return phi[1];
  }

  /* method for the drift rate */
  double drift(const double phi[9], double x, double t) const override {
    double mu = phi[2];
    double l = phi[3];
    double v = mu - l*x;
    return v;
  }

  /* method for the diffusion rate */
  double diffusion(const double phi[9], double x, double t) const override {
    return phi[4];
  }

  /* method for the upper threshold */
  double upper_threshold(const double phi[9], double t) const override {
    return phi[5];
  }

  /* method for the lower threshold */
  double lower_threshold(const double phi[9], double t) const override {
    return -phi[5];
  }

  /* method for the contamination strength */
  double contamination_strength(const double phi[9]) const override {
    return phi[6];
  }

  /* method for the contamination probability distribution */
  double contamination_probability(const double phi[9], double t) const override {
    double gl = phi[7];
    double gu = phi[8];
    double pg = 0.0;
    if ((t >= gl) && (t <= gu)) {
      pg = 1.0/(gu - gl);
    }
    return pg;
  }

  /* method for locally modifying the time step size */
  double modify_dt(const double phi[9], double t) const override {
    return 1.0;
  }

};



class UGM : public Model_TX {
protected:

  /* method for the non-decision time */
  double non_decision(const double phi[10]) const override {
    return phi[0];
  }

  /* method for the start point */
  double relative_start(const double phi[10]) const override {
    return phi[1];
  }

  /* method for the drift rate */
  double drift(const double phi[10], double x, double t) const override {
    double mu = phi[2];
    double l = pow(10.0, phi[3]);
    double k = pow(10.0, phi[4]);
    double v = mu*(1.0 + k*t) - (l - k/(1.0 + k*t) )*x;
    return v;
  }

  /* method for the diffusion rate */
  double diffusion(const double phi[10], double x, double t) const override {
    double k = pow(10.0, phi[4]);
    double sigma = phi[5];
    double D = sigma*(1.0 + k*t);
    return D;
  }

  /* method for the upper threshold */
  double upper_threshold(const double phi[10], double t) const override {
    return phi[6];
  }

  /* method for the lower threshold */
  double lower_threshold(const double phi[10], double t) const override {
    return -phi[6];
  }

  /* method for the contamination strength */
  double contamination_strength(const double phi[10]) const override {
    return phi[7];
  }

  /* method for the contamination probability distribution */
  double contamination_probability(const double phi[10], double t) const override {
    double gl = phi[8];
    double gu = phi[9];
    double pg = 0.0;
    if ((t >= gl) && (t <= gu)) {
      pg = 1.0/(gu - gl);
    }
    return pg;
  }

  /* method for locally modifying the time step size */
  double modify_dt(const double phi[10], double t) const override {
    return 1.0;
  }

};



class UGMF : public Model_TX {
protected:

  /* method for the non-decision time */
  double non_decision(const double phi[12]) const override {
    return phi[0];
  }

  /* method for the start point */
  double relative_start(const double phi[12]) const override {
    return phi[1];
  }

  /* method for the drift rate */
  double drift(const double phi[12], double x, double t) const override {
    double mu1 = phi[2];
    double mu2 = phi[3];
    double l = pow(10.0, phi[4]);
    double k = pow(10.0, phi[5]);
    double t0 = phi[6];
    double v = 0.0;

    if (t < t0) {
      v = mu1*(1.0 + k*t) - (l - k/(1.0 + k*t) )*x;
    } else {
      v = mu2*(1.0 + k*t) - (l - k/(1.0 + k*t) )*x;
    }

    return v;
  }

  /* method for the diffusion rate */
  double diffusion(const double phi[12], double x, double t) const override {
    double k = pow(10.0, phi[5]);
    double sigma = phi[7];
    double D = sigma*(1.0 + k*t);
    return D;
  }

  /* method for the upper threshold */
  double upper_threshold(const double phi[12], double t) const override {
    return phi[8];
  }

  /* method for the lower threshold */
  double lower_threshold(const double phi[12], double t) const override {
    return -phi[8];
  }

  /* method for the contamination strength */
  double contamination_strength(const double phi[12]) const override {
    return phi[9];
  }

  /* method for the contamination probability distribution */
  double contamination_probability(const double phi[12], double t) const override {
    double gl = phi[10];
    double gu = phi[11];
    double pg = 0.0;
    if ((t >= gl) && (t <= gu)) {
      pg = 1.0/(gu - gl);
    }
    return pg;
  }

  /* method for locally modifying the time step size */
  double modify_dt(const double phi[12], double t) const override {
    double t0 = phi[6];
    double range = 0.1;
    double dt_mod = 1.0;

    if ((t >= t0 - range) && (t <= t0 + range)) {
      dt_mod = 0.1;
    }

    return dt_mod;
  }

};


// ---- CUSTOM FUNCTION ----

// helper functions
static double callRFunction3(SEXP fun,
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
  double val = asReal(ans);
  UNPROTECT(3);
  return val;
}

static double callRFunction2(SEXP fun,
                             const double* phi, int n_phi, double t) {
  SEXP call, ans, phiR;
  PROTECT(phiR = Rf_allocVector(REALSXP, n_phi));
  for (int i = 0; i < n_phi; ++i)
    REAL(phiR)[i] = phi[i];
  PROTECT(call = Rf_lang3(fun, phiR, Rf_ScalarReal(t)));
  PROTECT(ans = Rf_eval(call, R_GlobalEnv));   // or another environment
  double val = asReal(ans);
  UNPROTECT(3);
  return val;
}

static double callRFunction1(SEXP fun, const double* phi, int n_phi) {
  SEXP call, ans, phiR;
  PROTECT(phiR = Rf_allocVector(REALSXP, n_phi));
  for (int i = 0; i < n_phi; ++i)
    REAL(phiR)[i] = phi[i];
  PROTECT(call = Rf_lang2(fun, phiR));
  PROTECT(ans = Rf_eval(call, R_GlobalEnv));   // or another environment
  double val = asReal(ans);
  UNPROTECT(3);
  return val;
}

// Structure holding optional user-defined functions for all overrideable methods
struct ModelTX_Callbacks {

  using Fn1   = double (*)(const double*, double, double);
  using Fn2   = double (*)(const double*, double);
  using Fn0   = double (*)(const double*);

  Fn1 drift                 = nullptr;
  Fn1 diffusion             = nullptr;
  Fn2 upper_threshold       = nullptr;
  Fn2 lower_threshold       = nullptr;
  Fn0 non_decision          = nullptr;
  Fn0 relative_start        = nullptr;
  Fn0 contamination_strength= nullptr;
  Fn2 contamination_probability = nullptr;
  Fn2 modify_dt             = nullptr;

  SEXP r_drift = R_NilValue;
  SEXP r_diffusion = R_NilValue;
  SEXP r_upper_threshold = R_NilValue;
  SEXP r_lower_threshold = R_NilValue;
  SEXP r_non_decision = R_NilValue;
  SEXP r_relative_start = R_NilValue;
  SEXP r_contamination_strength = R_NilValue;
  SEXP r_contamination_probability = R_NilValue;
  SEXP r_modify_dt = R_NilValue;

};

class CSTM_TX : public Model_TX {
public:
  static void set_callbacks(const ModelTX_Callbacks& cb) { callbacks = cb; }
  static const ModelTX_Callbacks& get_callbacks() { return callbacks; }

protected:

  /* method for the non-decision time */
  double non_decision(const double phi[100]) const override {
    if (callbacks.non_decision) {
      return callbacks.non_decision(phi);
    } else if (callbacks.r_non_decision) {
      return callRFunction1(callbacks.r_non_decision, phi, 100);
    } else {
      return phi[0];
    }
  }

  /* method for the start point */
  double relative_start(const double phi[100]) const override {
    if (callbacks.relative_start) {
      return callbacks.relative_start(phi);
    } else if (callbacks.r_relative_start) {
      return callRFunction1(callbacks.r_relative_start, phi, 100);
    } else {
      return phi[1];
    }
  }

  /* method for the drift rate */
  double drift(const double phi[100], double x, double t) const override {
    if (callbacks.drift) {
      return callbacks.drift(phi, x, t);
    } else if (callbacks.r_drift) {
      return callRFunction3(callbacks.r_drift, phi, 100, x, t);
    } else {
      return phi[2];
    }
  }

  /* method for the diffusion rate */
  double diffusion(const double phi[100], double x, double t) const override {
    if (callbacks.diffusion) {
      return callbacks.diffusion(phi, x, t);
    } else if (callbacks.r_diffusion) {
      return callRFunction3(callbacks.r_diffusion, phi, 100, x, t);
    } else {
      return phi[3];
    }
  }

  /* method for the upper threshold */
  double upper_threshold(const double phi[100], double t) const override {
    if (callbacks.upper_threshold) {
      return callbacks.upper_threshold(phi, t);
    } else if (callbacks.r_upper_threshold) {
      return callRFunction2(callbacks.r_upper_threshold, phi, 100, t);
    } else {
      return phi[4];
    }
  }

  /* method for the lower threshold */
  double lower_threshold(const double phi[100], double t) const override {
    if (callbacks.lower_threshold) {
      return callbacks.lower_threshold(phi, t);
    } else if (callbacks.r_lower_threshold) {
      return callRFunction2(callbacks.r_lower_threshold, phi, 100, t);
    } else {
      return -phi[4];
    }
  }

  /* method for the contamination strength */
  double contamination_strength(const double phi[100]) const override {
    if (callbacks.contamination_strength) {
      return callbacks.contamination_strength(phi);
    } else if (callbacks.r_contamination_strength) {
      return callRFunction1(callbacks.r_contamination_strength, phi, 100);
    } else {
      return phi[5];
    }
  }

  /* method for the contamination probability distribution */
  double contamination_probability(const double phi[100], double t) const override {
    if (callbacks.contamination_probability) {
      return callbacks.contamination_probability(phi, t);
    } else if (callbacks.r_contamination_probability) {
      return callRFunction2(callbacks.r_contamination_probability, phi, 100, t);
    } else {
      return (t >= phi[6] && t <= phi[7]) ? 1.0 / (phi[7] - phi[6]) : 0.0;
    }
  }

  /* method for locally modifying the time step size */
  double modify_dt(const double phi[100], double t) const override {
    if (callbacks.modify_dt) {
      return callbacks.modify_dt(phi, t);
    } else if (callbacks.r_modify_dt) {
      return callRFunction2(callbacks.r_modify_dt, phi, 100, t);
    } else {
      return 1.0;
    }
  }

private:
  static ModelTX_Callbacks callbacks;

};


#endif
