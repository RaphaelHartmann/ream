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



class CDSTP : public Model_TW {
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
  double drift_ts(const double phi[15]) const override {
    return phi[7];
  }

  /* function for the diffusion rate */
  double diffusion(const double phi[15], double t, double w) const override {
    double mu1 = phi[3] + phi[4]*phi[5];
    double mu2 = phi[6];
    double sigma = phi[8];
    double sigma_eff = phi[9];
    double D = sigma*sqrt(1.0 + sigma_eff*w);
    return D;
  }

  /* function for the target selection diffusion rate */
  double diffusion_ts(const double phi[15]) const override {
    return phi[8];
  }

  /* function for the upper threshold */
  double upper_threshold(const double phi[15], double t) const override {
    return phi[10];
  }

  /* function for the lower threshold */
  double lower_threshold(const double phi[15], double t) const override {
    return -phi[10];
  }

  /* function for the target selection upper threshold */
  double upper_threshold_ts(const double phi[15]) const override {
    return phi[11];
  }

  /* function for the target selection lower threshold */
  double lower_threshold_ts(const double phi[15]) const override {
    return -phi[11];
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
    double w_ts = relative_start_ts(phi); /* relative start point for process 2 */
    double v_ts = drift_ts(phi); /* drift rate for process 2 */
    double sigma_ts = diffusion_ts(phi); /* diffusion rate for process 2 */
    double a_ts = upper_threshold_ts(phi) - lower_threshold_ts(phi); /* threshold separation for process 2 */
    double z_ts = w_ts*a_ts; /* start point for process 2 */
    double v = 0.0; /* total summed drift rate for process 1 */
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
    double v = 0.0; /* total summed drift rate for process 1 */
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
    double mu1 = phi[3] + phi[4]*phi[5];
    double mu2 = phi[6];
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
    double mu1 = phi[3] + phi[4]*phi[5];
    double mu2 = phi[6];
    double lamb = phi[7];
    double kappa = phi[8];
    double w = 1.0 - exp( -pow(t/lamb, kappa) );
    return w;
  }
};



class CSTM_TW : public Model_TW {
protected:

  /* method for the non-decision time of process 1 */
  double non_decision(const double phi[100]) const override {
    return phi[0];
  }

  /* method for the start point of process 1 */
  double relative_start(const double phi[100]) const override {
    return phi[1];
  }

  /* method for the drift rate of process 1 */
  double drift(const double phi[100], double t, double w) const override {
    return phi[2];
  }

  /* method for the diffusion rate of process 1 */
  double diffusion(const double phi[100], double t, double w) const override {
    return phi[3];
  }

  /* method for the upper threshold of process 1 */
  double upper_threshold(const double phi[100], double t) const override {
    return phi[4];
  }

  /* method for the lower threshold of process 1 */
  double lower_threshold(const double phi[100], double t) const override {
    return -phi[4];
  }

  /* method for the contamination strength of process 1 */
  double contamination_strength(const double phi[100]) const override {
    return phi[5];
  }

  /* method for the contamination probability distribution of process 1 */
  double contamination_probability(const double phi[100], double t) const override {
    double gl = phi[6];
    double gu = phi[7];
    double pg = 0.0;
    if ((t >= gl) && (t <= gu)) {
      pg = 1.0/(gu - gl);
    }
    return pg;
  }

  /* method for locally modifying the time step size */
  double modify_dt(const double phi[100], double t) const override {
    return 1.0;
  }

};



#endif
