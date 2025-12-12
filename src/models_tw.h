/* file: model_functions.cpp
 Functions for defining model.
 Author: Mathew Murrow and Raphael Hartmann
 Date: Sep 02, 2024 */

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

#ifndef MODELS_TW_H
#define MODELS_TW_H

#include "Model_TW.h"
#include <Rinternals.h>


// class CDSTP : public Model_TW {
// protected:
//
//   /* function for the non-decision time */
//   double non_decision(const double phi[15]) const override {
//     return phi[0];
//   }
//
//   /* function for the start point */
//   double relative_start(const double phi[15]) const override {
//     return phi[1];
//   }
//
//   /* function for the target selection start point */
//   double relative_start_ts(const double phi[15]) const override {
//     return phi[2];
//   }
//
//   /* function for the drift rate */
//   double drift(const double phi[15], double t, double w) const override {
//     double mu1 = phi[3] + phi[4]*phi[5];
//     double mu2 = phi[6];
//     double v = (1.0 - w)*mu1 + w*mu2;
//     return v;
//   }
//
//   /* function for the target selection drift rate */
//   double drift_ts(const double phi[15]) const override {
//     return phi[7];
//   }
//
//   /* function for the diffusion rate */
//   double diffusion(const double phi[15], double t, double w) const override {
//     double sigma = phi[8];
//     double sigma_eff = phi[9];
//     double D = sigma*sqrt(1.0 + sigma_eff*w);
//     return D;
//   }
//
//   /* function for the target selection diffusion rate */
//   double diffusion_ts(const double phi[15]) const override {
//     return phi[8];
//   }
//
//   /* function for the upper threshold */
//   double upper_threshold(const double phi[15], double t) const override {
//     return phi[10];
//   }
//
//   /* function for the lower threshold */
//   double lower_threshold(const double phi[15], double t) const override {
//     return -phi[10];
//   }
//
//   /* function for the target selection upper threshold */
//   double upper_threshold_ts(const double phi[15]) const override {
//     return phi[11];
//   }
//
//   /* function for the target selection lower threshold */
//   double lower_threshold_ts(const double phi[15]) const override {
//     return -phi[11];
//   }
//
//   /* function for the contamination strength */
//   double contamination_strength(const double phi[15]) const override {
//     return phi[12];
//   }
//
//   /* function for the contamination probability distribution */
//   double contamination_probability(const double phi[15], double t) const override {
//     double gl = phi[13];
//     double gu = phi[14];
//     double pg = 0.0;
//     if ((t >= gl) && (t <= gu)) {
//       pg = 1.0/(gu - gl);
//     }
//     return pg;
//   }
//
//   /* function for locally modifying the time step size */
//   double modify_dt(const double phi[15], double t) const override {
//     return 1.0;
//   }
//
//   /* function used to calculate the CDF of the target selection process, used to set w(t) */
//   double ts_cdf(const double phi[15], double t) const override {
//     double w_ts = relative_start_ts(phi); /* relative start point for process 2 */
//     double v_ts = drift_ts(phi); /* drift rate for process 2 */
//     double sigma_ts = diffusion_ts(phi); /* diffusion rate for process 2 */
//     double a_ts = upper_threshold_ts(phi) - lower_threshold_ts(phi); /* threshold separation for process 2 */
//     double z_ts = w_ts*a_ts; /* start point for process 2 */
//     int kk = 0; /* looping index  */
//     int N_k = 0; /* number of iterations in infinite sum */
//
//     /* set number of iterations in infinite sum */
//     if (t <= flip) {
//       N_k = its_smalltime;
//     } else {
//       N_k = its_bigtime;
//     }
//
//     /* calculate probability p of process 2 crossing upper and lower threhsolds */
//     double p_lower = ( exp(-2.0*v_ts*a_ts/(sigma_ts*sigma_ts)) - exp(-2.0*v_ts*z_ts/(sigma_ts*sigma_ts)) ) / ( exp(-2.0*v_ts*a_ts/(sigma_ts*sigma_ts)) - 1.0 );
//
//     /* calculate cumulative probability g for upper and lower threhsolds */
//     double g_lower = 0.0;
//     for (kk = 1; kk < N_k; kk++) {
//       g_lower += 2.0*kk*sin(kk*pi*z_ts/a_ts)*exp(-0.5*t*((v_ts*v_ts)/(sigma_ts*sigma_ts) + (pi*pi)*(kk*kk)*(sigma_ts*sigma_ts)/(a_ts*a_ts))) / ((v_ts*v_ts)/(sigma_ts*sigma_ts) + (pi*pi)*(kk*kk)*(sigma_ts*sigma_ts)/(a_ts*a_ts));
//     }
//     g_lower = p_lower - pi*(sigma_ts*sigma_ts)/(a_ts*a_ts)*exp(-v_ts*z_ts/(sigma_ts*sigma_ts))*g_lower;
//
//     /* calculate w(t) */
//     double weight = g_lower/p_lower;
//     if (weight < 0.0) {
//       weight = 0.0;
//     }
//     if (weight > 1.0){
//       weight = 1.0;
//     }
//
//     return weight;
//   }
//
// };



class SDPM : public Model_TW {
protected:

  /* function for the non-decision time */
  double non_decision(const double phi[12]) const override {
    return phi[0];
  }

  /* function for the start point */
  double relative_start(const double phi[12]) const override {
    return phi[1]*phi[2];
  }

  /* function for the target selection start point */
  double relative_start_ts(const double phi[12]) const override {
    return phi[2];
  }

  /* function for the drift rate */
  double drift(const double phi[12], double t, double w) const override {
    double mu2 = phi[3];
    double v = w*mu2;
    return v;
  }

  /* function for the target selection drift rate */
  double drift_ts(const double phi[12]) const override {
    return phi[4];
  }

  /* function for the diffusion rate */
  double diffusion(const double phi[12], double t, double w) const override {
    double sigma = phi[5];
    double sigma_eff = phi[6];
    double D = sigma*sqrt(1.0 + sigma_eff*w);
    return D;
  }

  /* function for the target selection diffusion rate */
  double diffusion_ts(const double phi[12]) const override {
    return phi[5];
  }

  /* function for the upper threshold */
  double upper_threshold(const double phi[12], double t) const override {
    return phi[7];
  }

  /* function for the lower threshold */
  double lower_threshold(const double phi[12], double t) const override {
    return -phi[7];
  }

  /* function for the target selection upper threshold */
  double upper_threshold_ts(const double phi[12]) const override {
    return phi[8];
  }

  /* function for the target selection lower threshold */
  double lower_threshold_ts(const double phi[12]) const override {
    return -phi[8];
  }

  /* function for the contamination strength */
  double contamination_strength(const double phi[12]) const override {
    return phi[9];
  }

  /* function for the contamination probability distribution */
  double contamination_probability(const double phi[12], double t) const override {
    double gl = phi[10];
    double gu = phi[11];
    double pg = 0.0;
    if ((t >= gl) && (t <= gu)) {
      pg = 1.0/(gu - gl);
    }
    return pg;
  }

  /* function for locally modifying the time step size */
  double modify_dt(const double phi[12], double t) const override {
    return 1.0;
  }

  /* function used to calculate the CDF of the target selection process, used to set w(t) */
  double ts_cdf(const double phi[12], double t) const override {
    double w_ts = relative_start_ts(phi); /* relative start point for process 2 */
    double v_ts = drift_ts(phi); /* drift rate for process 2 */
    double sigma_ts = diffusion_ts(phi); /* diffusion rate for process 2 */
    double a_ts = upper_threshold_ts(phi) - lower_threshold_ts(phi); /* threshold separation for process 2 */
    double z_ts = w_ts*a_ts; /* start point for process 2 */
    int kk = 0; /* looping index  */
    int N_k = 0; /* number of iterations in infinite sum */

    /* set number of iterations in infinite sum */
    if (t <= flip) {
      N_k = its_smalltime;
    } else {
      N_k = its_bigtime;
    }

    /* calculate probability p of process 2 crossing upper and lower threhsolds */
    double p_lower = ( exp(-2.0*v_ts*a_ts/(sigma_ts*sigma_ts)) - exp(-2.0*v_ts*z_ts/(sigma_ts*sigma_ts)) ) / ( exp(-2.0*v_ts*a_ts/(sigma_ts*sigma_ts)) - 1.0 );

    /* calculate cumulative probability g for upper and lower threhsolds */
    double g_lower = 0.0;
    for (kk = 1; kk < N_k; kk++) {
      g_lower += 2.0*kk*sin(kk*pi*z_ts/a_ts)*exp(-0.5*t*((v_ts*v_ts)/(sigma_ts*sigma_ts) + (pi*pi)*(kk*kk)*(sigma_ts*sigma_ts)/(a_ts*a_ts))) / ((v_ts*v_ts)/(sigma_ts*sigma_ts) + (pi*pi)*(kk*kk)*(sigma_ts*sigma_ts)/(a_ts*a_ts));
    }
    g_lower = p_lower - pi*(sigma_ts*sigma_ts)/(a_ts*a_ts)*exp(-v_ts*z_ts/(sigma_ts*sigma_ts))*g_lower;

    /* calculate w(t) */
    double weight = g_lower/p_lower;
    if (weight < 0.0) {
      weight = 0.0;
    }
    if (weight > 1.0){
      weight = 1.0;
    }

    return weight;
  }

};



class WDSTP : public Model_TW {
protected:

  /* function for the non-decision time */
  double non_decision(const double phi[15]) const override {
    return phi[0];
  }

  /* function for the start point */
  double relative_start(const double phi[15]) const override {
    return phi[1];
  }

  /* function for the target selection start point */
  double relative_start_ts(const double phi[15]) const override {
    return phi[2];
  }

  /* function for the drift rate */
  double drift(const double phi[15], double t, double w) const override {
    double mu1 = phi[3] + phi[4]*phi[5];
    double mu2 = phi[6];
    double v = (1.0 - w)*mu1 + w*mu2;
    return v;
  }

  /* function for the target selection drift rate */
  double drift_ts(const double phi [15]) const override {
    return 0;
  }

  /* function for the diffusion rate */
  double diffusion(const double phi[15], double t, double w) const override {
    double sigma = phi[9];
    double sigma_eff = phi[10];
    double D = sigma*sqrt(1.0 + sigma_eff*w);
    return D;
  }

  /* function for the target selection diffusion rate */
  double diffusion_ts(const double phi[15]) const override {
    return 0;
  }

  /* function for the upper threshold */
  double upper_threshold(const double phi[15], double t) const override {
    return phi[11];
  }

  /* function for the lower threshold */
  double lower_threshold(const double phi[15], double t) const override {
    return -phi[11];
  }

  /* function for the target selection upper threshold */
  double upper_threshold_ts(const double phi[12]) const override {
    return 0;
  }

  /* function for the target selection lower threshold */
  double lower_threshold_ts(const double phi[12]) const override {
    return 0;
  }

  /* function for the contamination strength */
  double contamination_strength(const double phi[15]) const override {
    return phi[12];
  }

  /* function for the contamination probability distribution */
  double contamination_probability(const double phi[15], double t) const override {
    double gl = phi[13];
    double gu = phi[14];
    double pg = 0.0;
    if ((t >= gl) && (t <= gu)) {
      pg = 1.0/(gu - gl);
    }
    return pg;
  }

  /* function for locally modifying the time step size */
  double modify_dt(const double phi[15], double t) const override {
    return 1.0;
  }

  /* function used to calculate the CDF of the target selection process, used to set w(t) */
  double ts_cdf(const double phi[15], double t) const override {
    double lamb = phi[7];
    double kappa = phi[8];
    double w = 1.0 - exp( -pow(t/lamb, kappa) );
    return w;
  }
};


// ---- CUSTOM FUNCTION ----

// Structure holding optional user-defined functions for all overrideable methods
struct ModelTW_Callbacks {

  using Fn3   = double (*)(const double*, double, double);
  using Fn2   = double (*)(const double*, double);
  using Fn1   = double (*)(const double*);

  Fn3 drift                 = nullptr;
  Fn1 drift_ts              = nullptr;
  Fn3 diffusion             = nullptr;
  Fn1 diffusion_ts          = nullptr;
  Fn2 upper_threshold       = nullptr;
  Fn1 upper_threshold_ts    = nullptr;
  Fn2 lower_threshold       = nullptr;
  Fn1 lower_threshold_ts    = nullptr;
  Fn1 non_decision          = nullptr;
  Fn1 relative_start        = nullptr;
  Fn1 relative_start_ts     = nullptr;
  Fn1 contamination_strength= nullptr;
  Fn2 contamination_probability = nullptr;
  Fn2 modify_dt             = nullptr;
  Fn2 ts_cdf                = nullptr;

  SEXP r_drift = R_NilValue;
  SEXP r_drift_ts = R_NilValue;
  SEXP r_diffusion = R_NilValue;
  SEXP r_diffusion_ts = R_NilValue;
  SEXP r_upper_threshold = R_NilValue;
  SEXP r_upper_threshold_ts = R_NilValue;
  SEXP r_lower_threshold = R_NilValue;
  SEXP r_lower_threshold_ts = R_NilValue;
  SEXP r_non_decision = R_NilValue;
  SEXP r_relative_start = R_NilValue;
  SEXP r_relative_start_ts = R_NilValue;
  SEXP r_contamination_strength = R_NilValue;
  SEXP r_contamination_probability = R_NilValue;
  SEXP r_modify_dt = R_NilValue;
  SEXP r_ts_cdf = R_NilValue;

};

class CSTM_TW : public Model_TW {
public:
  static void set_callbacks(const ModelTW_Callbacks& cb) { callbacks = cb; }
  static ModelTW_Callbacks& get_callbacks() { return callbacks; }

protected:

  /* method for the non-decision time of process 1 */
  double non_decision(const double phi[100]) const override {
    if (callbacks.non_decision) {
      return callbacks.non_decision(phi);
    } else if (callbacks.r_non_decision != R_NilValue) {
      return callRFunction1(callbacks.r_non_decision, phi, 100);
    } else {
      return phi[0];
    }
  }

  /* method for the start point of process 1 */
  double relative_start(const double phi[100]) const override {
    if (callbacks.relative_start) {
      return callbacks.relative_start(phi);
    } else if (callbacks.r_relative_start != R_NilValue) {
      return callRFunction1(callbacks.r_relative_start, phi, 100);
    } else {
      return phi[1];
    }
  }

  /* function for the target selection start point */
  double relative_start_ts(const double phi[100]) const override {
    if (callbacks.relative_start_ts) {
      return callbacks.relative_start_ts(phi);
    } else if (callbacks.r_relative_start_ts != R_NilValue) {
      return callRFunction1(callbacks.r_relative_start_ts, phi, 100);
    } else {
      return 0.0;
    }
  }

  /* method for the drift rate of process 1 */
  double drift(const double phi[100], double t, double w) const override {
    if (callbacks.drift) {
      return callbacks.drift(phi, t, w);
    } else if (callbacks.r_drift != R_NilValue) {
      return callRFunction3(callbacks.r_drift, phi, 100, t, w);
    } else {
      return phi[2];
    }
  }

  /* function for the target selection drift rate */
  double drift_ts(const double phi[100]) const override {
    if (callbacks.drift_ts) {
      return callbacks.drift_ts(phi);
    } else if (callbacks.r_drift_ts != R_NilValue) {
      return callRFunction1(callbacks.r_drift_ts, phi, 100);
    } else {
      return 0.0;
    }
  }

  /* method for the diffusion rate of process 1 */
  double diffusion(const double phi[100], double t, double w) const override {
    if (callbacks.diffusion) {
      return callbacks.diffusion(phi, t, w);
    } else if (callbacks.r_diffusion != R_NilValue) {
      return callRFunction3(callbacks.r_diffusion, phi, 100, t, w);
    } else {
      return phi[3];
    }
  }

  /* function for the target selection diffusion rate */
  double diffusion_ts(const double phi[100]) const override {
    if (callbacks.diffusion_ts) {
      return callbacks.diffusion_ts(phi);
    } else if (callbacks.r_diffusion_ts != R_NilValue) {
      return callRFunction1(callbacks.r_diffusion_ts, phi, 100);
    } else {
      return phi[3];
    }
  }

  /* method for the upper threshold of process 1 */
  double upper_threshold(const double phi[100], double t) const override {
    if (callbacks.upper_threshold) {
      return callbacks.upper_threshold(phi, t);
    } else if (callbacks.r_upper_threshold != R_NilValue) {
      return callRFunction2(callbacks.r_upper_threshold, phi, 100, t);
    } else {
      return phi[4];
    }
  }

  /* method for the lower threshold of process 1 */
  double lower_threshold(const double phi[100], double t) const override {
    if (callbacks.lower_threshold) {
      return callbacks.lower_threshold(phi, t);
    } else if (callbacks.r_lower_threshold != R_NilValue) {
      return callRFunction2(callbacks.r_lower_threshold, phi, 100, t);
    } else {
      return -phi[4];
    }
  }

  /* function for the target selection upper threshold */
  double upper_threshold_ts(const double phi[100]) const override {
    if (callbacks.upper_threshold_ts) {
      return callbacks.upper_threshold_ts(phi);
    } else if (callbacks.r_upper_threshold_ts != R_NilValue) {
      return callRFunction1(callbacks.r_upper_threshold_ts, phi, 100);
    } else {
      return phi[4];
    }
  }

  /* function for the target selection lower threshold */
  double lower_threshold_ts(const double phi[100]) const override {
    if (callbacks.lower_threshold_ts) {
      return callbacks.lower_threshold_ts(phi);
    } else if (callbacks.r_lower_threshold_ts != R_NilValue) {
      return callRFunction1(callbacks.r_lower_threshold_ts, phi, 100);
    } else {
      return -phi[4];
    }
  }

  /* method for the contamination strength */
  double contamination_strength(const double phi[100]) const override {
    if (callbacks.contamination_strength) {
      return callbacks.contamination_strength(phi);
    } else if (callbacks.r_contamination_strength != R_NilValue) {
      return callRFunction1(callbacks.r_contamination_strength, phi, 100);
    } else {
      return phi[5];
    }
  }

  /* method for the contamination probability distribution */
  double contamination_probability(const double phi[100], double t) const override {
    if (callbacks.contamination_probability) {
      return callbacks.contamination_probability(phi, t);
    } else if (callbacks.r_contamination_probability != R_NilValue) {
      return callRFunction2(callbacks.r_contamination_probability, phi, 100, t);
    } else {
      return (t >= phi[6] && t <= phi[7]) ? 1.0 / (phi[7] - phi[6]) : 0.0;
    }
  }

  /* method for locally modifying the time step size */
  double modify_dt(const double phi[100], double t) const override {
    if (callbacks.modify_dt) {
      return callbacks.modify_dt(phi, t);
    } else if (callbacks.r_modify_dt != R_NilValue) {
      return callRFunction2(callbacks.r_modify_dt, phi, 100, t);
    } else {
      return 1.0;
    }
  }

  /* function used to calculate the CDF of the target selection process, used to set w(t) */
  double ts_cdf(const double phi[100], double t) const override {
    if (callbacks.ts_cdf) {
      return callbacks.ts_cdf(phi, t);
    } else if (callbacks.r_ts_cdf != R_NilValue) {
      return callRFunction2(callbacks.r_ts_cdf, phi, 100, t);
    } else {
      return 1.0;
    }
  }

private:
  static ModelTW_Callbacks callbacks;

};



#endif
