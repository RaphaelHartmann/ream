/* file: model_functions.cpp
 Functions for defining model.
 Author: Raphael Hartmann
 Date: Sep 02, 2024 */

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

#ifndef MODEL_H
#define MODEL_H

#include "tools.h"
#include <memory>
#include <vector>
#include <R.h>
#include <R_ext/Random.h>
#include <Rmath.h>


/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */


/* Base model class */
class Model {
public:

  /* Constructor */
  Model();
  /* Virtual destructor */
  virtual ~Model();

  /* method used to calculate model PDF function */
  virtual int pdf(double *RsumlogPDF, double *RPDFlow, double *RPDFupp, double *RlogPDFlow, double *RlogPDFupp, std::vector<double> rtl, std::vector<double> rtu, double *phi) const=0;

  /* method used to calculate model CDF function */
  virtual int cdf(double *RsumlogCDF, double *RCDFlow, double *RCDFupp, double *RlogCDFlow, double *RlogCDFupp, std::vector<double> rtl, std::vector<double> rtu, double *phi) const=0;

  /* method used to draw random samples */
  virtual int rand(double *Rrt, double *phi) const=0;

  /* method used to construct grid for PDF */
  virtual int grid_pdf(double *Rrt, double *Rpdf_u, double *Rpdf_l, double *phi) const=0;

};





#endif
