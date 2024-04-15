

#ifndef TOOLS_H
#define TOOLS_H

#include "Model.h"
#include <cmath>
#include <iterator>
#include <vector>
#include <R.h>

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

/* Global variables */
extern const char *OUTPUT;
extern const char *INPUT;
extern const char *PHI;

extern const char *RTL; // for ll
extern const char *RTU; // for ll

extern const char *OUTPUT2; // for l
extern const char *OUTPUT3; // for l

extern const char *ModelName;

extern int N; // for simulation
extern double dt; // for simulation

extern int N_deps; // for ll
extern double dt_scale; // for ll
extern double rt_max; // for ll
extern int N_rtl; // for ll
extern int N_rtu; // for ll

extern int N_phi;

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

/* Functions */
double unif_L();

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

/* constants used by likelihood generating function */
const int N_deps_max = 501; /* number of spatial mesh points */
const int N_dt = 1001; /* maximum number of time steps */
const int N_dt_ini = 25; /* number of small, initial time steps to eliminate oscillations */
const double dt_ini_scale = 0.01; /* uset to set initial time steps, lower value gives smaller time steps */
const int N_dt_scaleup = 100; /* number of time steps to scale up from initial time step to normal time step */
const double dt_max = 0.1; /* maximum time step allowed for likelihood solver */
const double delt = 6.04e-6; /* numerical derivative step size */
const int N_decmax = 9; /* max number of time step decreases for changing thresholds */
const float ds_ratio_cutoff = 0.02; /* sets the max threshold collapse ratio */
const float dt_mod_scale = 0.5; /* sets the time step change if exceeds ds_ratio_cutoff */
const float threshold_cutoff = 1.0e-4; /* sets the minimum threshold value */
const float p_fpt_min = 1.0e-5; /* sets the minimum likelihood probability */
const float int_prob_min = 0.25; /* used to check if enough probability has accumulated to cutoff solver */

/* constants used by function approx_dt */
const double t_max = 100.0; /* simulate until this time only */
const int N_sims = 10; /* number of simulations */
const double dt_sims = 0.025; /* time step for simulation */

/* constants used by function simulate */
const double tsim_max = 100.0; /* simulate until this time only */

#endif
