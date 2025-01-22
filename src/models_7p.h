/* file: model_functions.cpp
 Functions for defining model.
 Author: Mathew Murrow and Raphael Hartmann
 Date: Jan 13, 2025 */

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

#ifndef MODELS_7P_H
#define MODELS_7P_H

#include "Model_7P.h"



class DDM : public Model_7P {
protected:

  /* function for the non-decision time's average */
  double non_decision_average(const double phi[11]) const override {
    return phi[0];
  }

  /* function for the non-decision time's width */
  double non_decision_width(const double phi[11]) const override {
    return phi[1];
  }

  /* function for the relative start's average */
  double relative_start_average(const double phi[11]) const override {
    return phi[2];
  }

  /* function for the relative start's width */
  double relative_start_width(const double phi[11]) const override {
    return phi[3];
  }

  /* function for the drift rate's average */
  double drift_average(const double phi[11]) const override {
    return phi[4];
  }

  /* function for the drift rate's width */
  double drift_width(const double phi[11]) const override {
    return phi[5];
  }

  /* function for the diffusion rate */
  double diffusion(const double phi[11], double x, double t) const override {
    return phi[6];
  }

  /* function for the upper threshold */
  double upper_threshold(const double phi[11], double t) const override {
    return phi[7];
  }

  /* function for the lower threshold */
  double lower_threshold(const double phi[11], double t) const override {
    return -phi[7];
  }

  /* function for the contamination strength */
  double contamination_strength(const double phi[11]) const override {
    return phi[8];
  }

  /* function for the contamination probability distribution */
  double contamination_probability(const double phi[11], double t) const override {
    double gl = phi[9];
    double gu = phi[10];
    double pg = 0.0;
    if ((t >= gl) && (t <= gu)) {
      pg = 1.0/(gu - gl);
    }
    return pg;
  }

  /* function for locally modifying the time step size */
  double modify_dt(const double phi[11], double t) const override {
    return 1.0;
  }

};



#endif
