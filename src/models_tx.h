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
    double l = pow(10.0, phi[4]);
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
    double l = pow(10.0, phi[3]);
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



class CSTM_TX : public Model_TX {
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
  double drift(const double phi[100], double x, double t) const override {
    return phi[2];
  }

  /* method for the diffusion rate of process 1 */
  double diffusion(const double phi[100], double x, double t) const override {
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
