/* file: model_functions.cpp
 Functions for defining model.
 Author: Raphael Hartmann and Mathew Murrow
 Date: Jan 13, 2025 */

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

#ifndef MODEL_7P_H
#define MODEL_7P_H

#include "Model.h"


class Model_7P : public Model {
public:
  /* Constructor */
  Model_7P();
  /* Virtual destructor */
  virtual ~Model_7P();

  /* method used to calculate model PDF function */
  int pdf(double *RsumlogPDF, double *RPDFlow, double *RPDFupp, double *RlogPDFlow, double *RlogPDFupp, std::vector<double> rtl, std::vector<double> rtu, double *phi) const override;

  /* method used to calculate model CDF function */
  int cdf(double *RsumlogCDF, double *RCDFlow, double *RCDFupp, double *RlogCDFlow, double *RlogCDFupp, std::vector<double> rtl, std::vector<double> rtu, double *phi) const override;

  /* method used to draw random samples */
  int rand(double *Rrt, double *phi) const override;

  /* method used to construct grid for PDF */
  int grid_pdf(double *Rrt, double *Rpdf_u, double *Rpdf_l, double *phi) const override;


protected:

  /* function for the non-decision time's average */
  virtual double non_decision_average(const double phi[]) const=0;

  /* function for the non-decision time's width */
  virtual double non_decision_width(const double phi[]) const=0;

  /* function for the relative start's average */
  virtual double relative_start_average(const double phi[]) const=0;

  /* function for the relative start's width */
  virtual double relative_start_width(const double phi[]) const=0;

  /* function for the drift rate's average */
  virtual double drift_average(const double phi[]) const=0;

  /* function for the drift rate's width */
  virtual double drift_width(const double phi[]) const=0;

  /* method for the diffusion rate */
  virtual double diffusion(const double phi[], double x, double t) const=0;

  /* method for the upper threshold */
  virtual double upper_threshold(const double phi[], double t) const=0;

  /* method for the lower threshold */
  virtual double lower_threshold(const double phi[], double t) const=0;

  /* method for the contamination strength */
  virtual double contamination_strength(const double phi[]) const=0;

  /* method for the contamination probability distribution */
  virtual double contamination_probability(const double phi[], double t) const=0;

  /* method for locally modifying the time step size */
  virtual double modify_dt(const double phi[], double t) const=0;

private:

  /* method which sets time step for likelihood generating function */
  double approx_dt(double* phi, double dt_scale) const;

  /* ---------------------------------------------------- */
  /* functions for calculating across-trial variabilities */
  /* ---------------------------------------------------- */

  /* pdf for standard normal distribution */
  double snd_pdf(double x) const;

  /* cdf for standard normal distribution */
  double snd_cdf(double x) const;

  /* function for interpolating an array */
  double interp(double t, double t0, double t1, double y0, double y1) const;

  /* uniform non-decision time distribution */
  double tnd_uni(double t, double a, double b) const;

  /* truncated normal non-decision time distribution */
  double tnd_tn(double* phi, double t) const;

  /* inverse gaussian (Wald) non-decision time distribution */
  double tnd_ig(double* phi, double t) const;

  /* uniform relative start point distribution */
  double w_uni(double eps, double a, double b) const;

  /* truncated normal relative start point distribution */
  double w_tn(double* phi, double eps, double a, double b) const;

  /* function for uniform drift rate distribution */
  double v_uni(double* phi, double v) const;

  /* function for normal drift rate distribution */
  double v_norm(double* phi, double v) const;

};

#endif
