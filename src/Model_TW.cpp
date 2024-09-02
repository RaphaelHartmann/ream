/* file: model_functions.cpp
 Functions for defining model.
 Author: Mathew Murrow and Raphael Hartmann
 Date: Sep 02, 2024 */

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

#include "Model_TW.h"


/* Model_TW Constructor */
Model_TW::Model_TW() {}

/* Model_TW Destructor */
Model_TW::~Model_TW() {}



/* method which sets time step for likelihood generating function */
double Model_TW::approx_dt(double* phi, double dt_scale) const {

  /* initialize variables */
  // int seed = 90210; /* set random number generator seed */
  int ii = 0, jj = 0; /* initialize loop indicies */
  double ww = relative_start(phi); /* relative start point */
  double vv = 0.0; /* drift rate */
  double DD = diffusion(phi, 0.0, 0.0); /* diffusion rate */
  double bu = upper_threshold(phi, 0.0); /* upper decision threshold */
  double bl = lower_threshold(phi, 0.0); /* lower decision threshold */
  double zz = bl + ww*(bu - bl); /* calculate absolute start point */
  double sqrtdt = sqrt(dt_sims); /* calculate square root of time step */
  double tt = 0.0; /* time */
  double xx = 0.0; /* accumulated evidence */
  double uu = 0.0; /* random number value */
  double fpt = 0.0; /* first passage time for threshold crossing */
  double dt_ = 0.0; /* time step output to likelihood generation function */
  double weight = 0.5;

  /* seed rng */
  // srand(seed);
  GetRNGstate();

  /* simulate model to determine time step */
  for (ii = 0; ii < N_sims; ii++) {

    /* reset time and space for next simulation */
    tt = 0.0;
    xx = zz;
    jj = 0;

    while ((jj < 1) && (tt <= t_max)) {

      /* update time step */
      tt += dt_sims;

      /* update drift rate, diffusion rate, and decision threhsolds */
      vv = drift(phi, tt, weight);
      DD = diffusion(phi, tt, weight);
      bu = upper_threshold(phi, tt);
      bl = lower_threshold(phi, tt);

      /* draw random number for simulating diffusion */
      // uu = -1 +  2 * ( rand() % 2 );
      uu = unif_L() >= 0.5 ? 1.0 : -1.0;

      /* calculate accumulated evidence */
      xx += dt_*vv + sqrtdt*DD*uu;

      /* check if accumulated evidence has crossed a decision threshold */
      if ((xx >= bu) || (xx <= bl)) {
        fpt += tt;
        jj = 1;
      }

    }
  }

  /* if no threshold crossing, set time step to max */
  if (fpt == 0.0) {
    fpt = N_sims*t_max;
  }

  /* calculate estimated time step for likelihood generation function */
  dt_ = dt_scale*fpt/N_sims;

  /* output time step for likelihood generation function */
  return dt_;

}


/* method for PDF */
int Model_TW::pdf(double *RsumlogPDF, double *RPDFlow, double *RPDFupp, double *RlogPDFlow, double *RlogPDFupp, std::vector<double> rtl, std::vector<double> rtu, double *phi) const {

  /* initialize loop indicies */
  int ii = 0, jj = 0, kk = 0;

  /* initialize arrays */
  double p_fpt[3][N_dt]; /* first passage time probability */
  double p_n[N_deps_max]; /* probability of evidence x at time step n */
  double p_ncn[N_deps_max]; /* Crank-Nicolson modified probability of evidence x at time step n */
  double p_np1[N_deps_max]; /* Crank-Nicolson modified probability of evidence x at time step n */
  double eps[N_deps_max]; /* probability of evidence x at time step n plus 1 */
  double AA[N_deps_max], BB[N_deps_max], CC[N_deps_max], DD[N_deps_max], EE[N_deps_max], FF[N_deps_max]; /* Crank-Nicolson matrix arrays */
  double ws[2]; /* start point array */
  int wi[2]; /* index of start point arrays */

  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < N_dt; jj++) {
      p_fpt[ii][jj] = 0.0;
    }
  }

  for (ii = 0; ii < N_deps_max; ii++) {
    p_n[ii] = 0.0;
    p_ncn[ii] = 0.0;
    p_np1[ii] = 0.0;
    eps[ii] = 0.0;
    AA[ii] = 0.0;
    BB[ii] = 0.0;
    CC[ii] = 0.0;
    DD[ii] = 0.0;
    EE[ii] = 0.0;
    FF[ii] = 0.0;
  }

  wi[0] = 0.0;
  wi[1] = 0.0;
  ws[0] = 0.0;
  ws[1] = 0.0;

  /* initilize constants */
  double v_n = 0.0, v_np1 = 0.0; /* drift rate at time steps n and n plus 1 */
  double sigma_n = 0.0, sigma_np1 = 0.0;  /* diffusion rate at time steps n and n plus 1 */
  double sigma2_n = 0.0, sigma2_np1 = 0.0; /* squared diffusion rate at time steps n and n plus 1 */
  double deps = 1.0/(N_deps - 1.0);
  double dt_base = 0.0, dt_ini = 0.0, dt_ = 0.0;
  double t_in = 0.0;
  double bu_nm1 = 0.0, bu_n = 0.0, bu_np1 = 0.0;
  double bl_nm1 = 0.0, bl_n = 0.0, bl_np1 = 0.0;
  double s_nm1 = 0.0, s_n = 0.0, s_np1 = 0.0;
  double dbudt_n = 0.0, dbudt_np1 = 0.0;
  double dbldt_n = 0.0, dbldt_np1 = 0.0;
  double dsdt_n = 0.0, dsdt_np1 = 0.0;
  double tnd = non_decision(phi);
  double ww = relative_start(phi);
  double tt = 0.0;
  double ds_ratio = 0.0, ds_ratio_np1 = 0.0, ds_ratio_nm1 = 0.0;
  double x_n = 0.0, x_np1 = 0.0;
  double alpha_n = 0.0, alpha_np1 = 0.0;
  double beta_n = 0.0, beta_np1 = 0.0;
  double gamma_n = 0.0, gamma_np1 = 0.0;
  double WW = 0.0;
  double s_ini = upper_threshold(phi, 0.0) - lower_threshold(phi, 0.0);
  double int_prob = 0.0;
  double gg = contamination_strength(phi);
  double g_prob = 0.0;
  double weight = 0.0;
  int N_cut = N_dt;

  /* set value of first element for first passage time array */
  p_fpt[0][0] = tnd;
  g_prob = 0.5*gg*contamination_probability(phi, p_fpt[0][0]);
  p_fpt[1][0] = g_prob;
  p_fpt[2][0] = g_prob;

  /* create array of spatial locations */
  eps[0] = 0.0;
  for (ii = 0; ii < N_deps; ii++) {
    eps[ii] = ii*deps;
  }
  eps[N_deps-1] = 1.0;

  /* calculate time step sizes */
  dt_base = approx_dt(phi, dt_scale);
  if (dt_base > dt_max) {
    dt_base = dt_max;
  }
  dt_ini = dt_ini_scale*dt_base;
  dt_ = dt_ini;

  /* calculate initial threshold locations */
  t_in = dt_ini/10.0;
  bu_np1 = upper_threshold(phi, t_in);
  bl_np1 = lower_threshold(phi, t_in);

  /* calculate initial threshold derivatives */
  dbudt_np1 = ( upper_threshold(phi, t_in+delt) - upper_threshold(phi, t_in-delt) ) / (2.0*delt);
  dbldt_np1 = ( lower_threshold(phi, t_in+delt) - lower_threshold(phi, t_in-delt) ) / (2.0*delt);

  /* calculate initial threhsold separation and derivative */
  s_np1 = bu_np1 - bl_np1;
  dsdt_np1 = dbudt_np1 - dbldt_np1;

  /* calculate probability around start point */
  /* determine nearest array elements to delta distribution */
  wi[0] = int( ceil(ww/deps) );
  wi[1] = int( floor(ww/deps) );

  /* determine distance between delta distribution and nearest array elements */
  ws[0] = fabs( ww - deps*wi[0] );
  ws[1] = fabs( ww - deps*wi[1] );

  /* numerically approximate delta function, distribute between nearest array elements */
  if (wi[0] == wi[1]) {
    p_n[ wi[0] ] = 1.0/deps;
  } else {
    p_n[ wi[0] ] = (1.0 - ws[0]/deps)/deps;
    p_n[ wi[1] ] = (1.0 - ws[1]/deps)/deps;
  }

  /* calculate first passage time distribution */
  for (ii = 0; ii < N_dt - 1; ii++) {

    /* set decision threshold location at previous time step*/
    bu_n = bu_np1;
    dbudt_n = dbudt_np1;
    bl_n = bl_np1;
    dbldt_n = dbldt_np1;
    s_n = bu_n - bl_n;
    dsdt_n = dbudt_n - dbldt_n;

    /* set time step */
    if (ii < N_dt_ini) {
      dt_ = dt_ini;
    } else if ((ii >= N_dt_ini) && (ii < N_dt_ini + N_dt_scaleup)) {
      dt_ = dt_ini + (ii + 1.0 - N_dt_ini)*(dt_base - dt_ini)/(N_dt_scaleup);
    } else {
      dt_ = dt_base;
    }

    /* decrease time step if decision thresholds decrease too rapidly */
    jj = 0;
    kk = 0;
    while ((jj < 1) && (kk < N_decmax)) {
      t_in = tt + dt_;
      bu_np1 = upper_threshold(phi, t_in);
      bl_np1 = lower_threshold(phi, t_in);

      t_in = tt - dt_;
      if (t_in <= 0.0) {
        t_in = dt_ini/10.0;
        bu_nm1 = upper_threshold(phi, t_in);
        bl_nm1 = lower_threshold(phi, t_in);
      } else {
        bu_nm1 = upper_threshold(phi, t_in);
        bl_nm1 = lower_threshold(phi, t_in);
      }

      s_np1 = bu_np1 - bl_np1;
      s_nm1 = bu_nm1 - bl_nm1;

      ds_ratio = 0.0;
      ds_ratio_np1 = fabs(s_n - s_np1)/s_n;
      ds_ratio_nm1 = fabs(s_n - s_nm1)/s_n;

      if (ds_ratio_np1 < ds_ratio_nm1) {
        ds_ratio = ds_ratio_nm1;
      } else {
        ds_ratio = ds_ratio_np1;
      }

      if (ds_ratio > ds_ratio_cutoff) {
        dt_ = dt_*dt_mod_scale;
      } else {
        jj = 1;
      }

      kk += 1;

    }

    /* calculate new decision thresholds */
    dt_ = modify_dt(phi, tt)*dt_;
    tt += dt_;

    bu_np1 = upper_threshold(phi, tt);
    dbudt_np1 = ( upper_threshold(phi, tt+delt) - upper_threshold(phi, tt-delt) ) / (2.0*delt);

    if (bu_np1 <= threshold_cutoff) {
      bu_np1 = threshold_cutoff;
      dbudt_np1 = 0.0;
    }

    bl_np1 = lower_threshold(phi, tt);
    dbldt_np1 = ( lower_threshold(phi, tt+delt) - lower_threshold(phi, tt-delt) ) / (2.0*delt);

    if (bl_np1 >= -threshold_cutoff) {
      bl_np1 = -threshold_cutoff;
      dbldt_np1 = 0.0;
    }

    s_np1 = bu_np1 - bl_np1;
    dsdt_np1 = dbudt_np1 - dbldt_np1;

    /* calculate drift and diffusion rates */
    weight = ts_cdf(phi, tt);

    if (ii == 0) {
      v_n = drift(phi, tt-dt_, weight);
      sigma_n = diffusion(phi, tt-dt_, weight);
      sigma2_n = sigma_n*sigma_n;
    } else {
      sigma_n = sigma_np1;
      sigma2_n = sigma2_np1;
      v_n = v_np1;
    }
    v_np1 = drift(phi, tt, weight);
    sigma_np1 = diffusion(phi, tt, weight);
    sigma2_np1 = sigma_np1*sigma_np1;

    /* invert matrix using tridiagonal matrix algorithm */
    for (jj = 0; jj < 2; jj++) {
      p_np1[jj] = 0.0;
    }

    // for (jj = 1; jj < 2; jj++) {
    //
    //   p_np1[jj] = 0.0;
    //
    //   if (ii == 0) {
    //     x_n = s_n*eps[jj] + bl_n;
    //     sigma_n[jj] = diffusion(phi, x_n, tt-dt_);
    //     sigma2_n[jj] = sigma_n[jj]*sigma_n[jj];
    //   } else {
    //     x_n = x_np1;
    //     sigma_n[jj] = sigma_np1[jj];
    //     sigma2_n[jj] = sigma2_np1[jj];
    //   }
    //
    //   x_np1 = s_np1*eps[jj] + bl_np1;
    //   sigma_np1[jj] = diffusion(phi, x_np1, tt);
    //   sigma2_np1[jj] = sigma_np1[jj]*sigma_np1[jj];
    //
    // }

    for (jj = 1; jj < N_deps-1; jj++) {

      p_np1[jj+1] = 0.0;

      alpha_n = dt_/(4.0*s_n*deps);
      alpha_np1 = dt_/(4.0*s_np1*deps);

      beta_n = eps[jj]*dsdt_n + dbldt_n;
      beta_np1 = eps[jj]*dsdt_np1 + dbldt_np1;

      gamma_n = alpha_n/(s_n*deps);
      gamma_np1 = alpha_np1/(s_np1*deps);

      AA[jj] = alpha_np1*beta_np1 - alpha_np1*v_np1 - gamma_np1*sigma2_np1;
      BB[jj] = 1.0 + 2.0*gamma_np1*sigma2_np1;
      CC[jj] = -alpha_np1*beta_np1 + alpha_np1*v_np1 - gamma_np1*sigma2_np1;

      DD[jj] = -alpha_n*beta_n + alpha_n*v_n + gamma_n*sigma2_n;
      EE[jj] = 1.0 - 2.0*gamma_n*sigma2_n;
      FF[jj] = alpha_n*beta_n - alpha_n*v_n + gamma_n*sigma2_n;

      p_ncn[jj] = DD[jj]*p_n[jj-1] + EE[jj]*p_n[jj] + FF[jj]*p_n[jj+1];

    }

    for (jj = 2; jj < N_deps - 1; jj++) {
      WW = AA[jj]/BB[jj-1];
      BB[jj] = BB[jj] - WW*CC[jj-1];
      p_ncn[jj] = p_ncn[jj] - WW*p_ncn[jj-1];
    }

    p_np1[N_deps-2] = p_ncn[N_deps-2]/BB[N_deps-2];
    p_n[N_deps-2] = p_np1[N_deps-2];

    for (jj = N_deps - 3; jj > 0; jj--) {
      p_np1[jj] = (p_ncn[jj] - CC[jj]*p_np1[jj+1])/BB[jj];
      p_n[jj] = p_np1[jj];
    }

    /* calculate probability of first passage time */
    p_fpt[0][ii+1] = tt + tnd;
    p_fpt[1][ii+1] = fabs( 0.5 * sigma2_np1 * (4.0*p_np1[N_deps-2] - p_np1[N_deps-3])/(2.0*deps*s_np1*s_ini) );
    p_fpt[2][ii+1] = fabs( 0.5 * sigma2_np1 * (4.0*p_np1[1] - p_np1[2])/(2.0*deps*s_np1*s_ini) );

    /* if p_fpt less than p_fpt_min, set to p_fpt_min */
    if (p_fpt[1][ii+1] <= p_fpt_min) {
      p_fpt[1][ii+1] = p_fpt_min;
    }
    if (p_fpt[2][ii+1] <= p_fpt_min) {
      p_fpt[2][ii+1] = p_fpt_min;
    }

    /* determine integrated probability */
    int_prob += p_fpt[1][ii+1]*dt_ + p_fpt[2][ii+1]*dt_;

    /* stop solver if time above rt_max or probability below p_fpt_min */
    if (p_fpt[0][ii+1] > rt_max) {
      g_prob = 0.5*gg*contamination_probability(phi, p_fpt[0][ii+1]);
      p_fpt[1][ii+1] = g_prob + (1.0 - gg)*p_fpt[1][ii+1];
      p_fpt[2][ii+1] = g_prob + (1.0 - gg)*p_fpt[2][ii+1];
      N_cut = ii + 1;
      break;
    } else if ((int_prob > int_prob_min) && (p_fpt[1][ii+1] <= p_fpt_min) && (p_fpt[2][ii+1] <= p_fpt_min)) {
      g_prob = 0.5*gg*contamination_probability(phi, p_fpt[0][ii+1]);
      p_fpt[1][ii+1] = g_prob + (1.0 - gg)*p_fpt[1][ii+1];
      p_fpt[2][ii+1] = g_prob + (1.0 - gg)*p_fpt[2][ii+1];
      N_cut = ii + 1;
      break;
    }

    /* calculate pdf including contamination distribution */
    g_prob = 0.5*gg*contamination_probability(phi, p_fpt[0][ii+1]);
    p_fpt[1][ii+1] = g_prob + (1.0 - gg)*p_fpt[1][ii+1];
    p_fpt[2][ii+1] = g_prob + (1.0 - gg)*p_fpt[2][ii+1];

  }

  /* determine size of rt vectors */
  int N_rtu = rtu.size();
  int N_rtl = rtl.size();

  /* initialize logPDF variable */
  double logPDF = 0.0;

  /* interpolate to find loglikelihood of data */
  kk = 0;
  for (ii = 0; ii < N_rtu; ii++){
    if ( (rtu[ii] >= p_fpt[0][N_cut]) || (rtu[ii] <= p_fpt[0][0]) ) {
      RPDFupp[ii] = p_fpt_min;
      RlogPDFupp[ii] = log(RPDFupp[ii]);
      logPDF += RlogPDFupp[ii];
    } else {
      for (jj = kk; jj < N_cut; jj++){
        if ( (rtu[ii] >= p_fpt[0][jj]) && (rtu[ii] < p_fpt[0][jj+1]) ) {
          RPDFupp[ii] = p_fpt[1][jj] + (rtu[ii] - p_fpt[0][jj]) * (p_fpt[1][jj+1] - p_fpt[1][jj]) / (p_fpt[0][jj+1] - p_fpt[0][jj]);
          RlogPDFupp[ii] = log(RPDFupp[ii]);
          logPDF += RlogPDFupp[ii];
          kk = jj;
          break;
        }
      }
    }
  }

  kk = 0;
  for (ii = 0; ii < N_rtl; ii++){
    if ( (rtl[ii] >= p_fpt[0][N_cut]) || (rtl[ii] <= p_fpt[0][0]) ) {
      RPDFlow[ii] = p_fpt_min;
      RlogPDFlow[ii] = log(RPDFlow[ii]);
      logPDF += RlogPDFlow[ii];
    } else {
      for (jj = kk; jj < N_cut; jj++){
        if ( (rtl[ii] >= p_fpt[0][jj]) && (rtl[ii] < p_fpt[0][jj+1]) ) {
          RPDFlow[ii] = p_fpt[2][jj] + (rtl[ii] - p_fpt[0][jj]) * (p_fpt[2][jj+1] - p_fpt[2][jj]) / (p_fpt[0][jj+1] - p_fpt[0][jj]);
          RlogPDFlow[ii] = log(RPDFlow[ii]);
          logPDF += RlogPDFlow[ii];
          kk = jj;
          break;
        }
      }
    }
  }

  /* write logPDF to R object */
  RsumlogPDF[0] = logPDF;

  return 0;

}


/* method used to calculate model CDF function */
int Model_TW::cdf(double *RsumlogCDF, double *RCDFlow, double *RCDFupp, double *RlogCDFlow, double *RlogCDFupp, std::vector<double> rtl, std::vector<double> rtu, double *phi) const {

  /* initialize loop indicies */
  int ii = 0, jj = 0, kk = 0;

  /* initialize arrays */
  double p_fpt[3][N_dt]; /* first passage time probability */
  double cdf[2][N_dt]; /* cumulative distribution function of first passage time probability */
  double p_n[N_deps_max]; /* probability of evidence x at time step n */
  double p_ncn[N_deps_max]; /* Crank-Nicolson modified probability of evidence x at time step n */
  double p_np1[N_deps_max]; /* Crank-Nicolson modified probability of evidence x at time step n */
  double eps[N_deps_max]; /* probability of evidence x at time step n plus 1 */
  double AA[N_deps_max], BB[N_deps_max], CC[N_deps_max], DD[N_deps_max], EE[N_deps_max], FF[N_deps_max]; /* Crank-Nicolson matrix arrays */
  double ws[2]; /* start point array */
  int wi[2]; /* index of start point arrays */

  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < N_dt; jj++) {
      p_fpt[ii][jj] = 0.0;
    }
  }

  for (ii = 0; ii < 2; ii++) {
    for (jj = 0; jj < N_dt; jj++) {
      cdf[ii][jj] = 0.0;
    }
  }

  for (ii = 0; ii < N_deps_max; ii++) {
    p_n[ii] = 0.0;
    p_ncn[ii] = 0.0;
    p_np1[ii] = 0.0;
    eps[ii] = 0.0;
    AA[ii] = 0.0;
    BB[ii] = 0.0;
    CC[ii] = 0.0;
    DD[ii] = 0.0;
    EE[ii] = 0.0;
    FF[ii] = 0.0;
  }

  wi[0] = 0.0;
  wi[1] = 0.0;
  ws[0] = 0.0;
  ws[1] = 0.0;

  /* initilize constants */
  double v_n = 0.0, v_np1 = 0.0; /* drift rate at time steps n and n plus 1 */
  double sigma_n = 0.0, sigma_np1 = 0.0;  /* diffusion rate at time steps n and n plus 1 */
  double sigma2_n = 0.0, sigma2_np1 = 0.0; /* squared diffusion rate at time steps n and n plus 1 */
  double deps = 1.0/(N_deps - 1.0);
  double dt_base = 0.0, dt_ini = 0.0, dt_ = 0.0;
  double t_in = 0.0;
  double bu_nm1 = 0.0, bu_n = 0.0, bu_np1 = 0.0;
  double bl_nm1 = 0.0, bl_n = 0.0, bl_np1 = 0.0;
  double s_nm1 = 0.0, s_n = 0.0, s_np1 = 0.0;
  double dbudt_n = 0.0, dbudt_np1 = 0.0;
  double dbldt_n = 0.0, dbldt_np1 = 0.0;
  double dsdt_n = 0.0, dsdt_np1 = 0.0;
  double tnd = non_decision(phi);
  double ww = relative_start(phi);
  double tt = 0.0;
  double ds_ratio = 0.0, ds_ratio_np1 = 0.0, ds_ratio_nm1 = 0.0;
  double x_n = 0.0, x_np1 = 0.0;
  double alpha_n = 0.0, alpha_np1 = 0.0;
  double beta_n = 0.0, beta_np1 = 0.0;
  double gamma_n = 0.0, gamma_np1 = 0.0;
  double WW = 0.0;
  double s_ini = upper_threshold(phi, 0.0) - lower_threshold(phi, 0.0);
  double int_prob = 0.0;
  double gg = contamination_strength(phi);
  double g_prob = 0.0;
  double weight = 0.0;
  int N_cut = N_dt;

  /* set value of first element for first passage time array */
  p_fpt[0][0] = tnd;
  g_prob = 0.5*gg*contamination_probability(phi, p_fpt[0][0]);
  p_fpt[1][0] = g_prob;
  p_fpt[2][0] = g_prob;

  /* create array of spatial locations */
  eps[0] = 0.0;
  for (ii = 0; ii < N_deps; ii++) {
    eps[ii] = ii*deps;
  }
  eps[N_deps-1] = 1.0;

  /* calculate time step sizes */
  dt_base = approx_dt(phi, dt_scale);
  if (dt_base > dt_max) {
    dt_base = dt_max;
  }
  dt_ini = dt_ini_scale*dt_base;
  dt_ = dt_ini;

  /* calculate initial threshold locations */
  t_in = dt_ini/10.0;
  bu_np1 = upper_threshold(phi, t_in);
  bl_np1 = lower_threshold(phi, t_in);

  /* calculate initial threshold derivatives */
  dbudt_np1 = ( upper_threshold(phi, t_in+delt) - upper_threshold(phi, t_in-delt) ) / (2.0*delt);
  dbldt_np1 = ( lower_threshold(phi, t_in+delt) - lower_threshold(phi, t_in-delt) ) / (2.0*delt);

  /* calculate initial threhsold separation and derivative */
  s_np1 = bu_np1 - bl_np1;
  dsdt_np1 = dbudt_np1 - dbldt_np1;

  /* calculate probability around start point */
  /* determine nearest array elements to delta distribution */
  wi[0] = int( ceil(ww/deps) );
  wi[1] = int( floor(ww/deps) );

  /* determine distance between delta distribution and nearest array elements */
  ws[0] = fabs( ww - deps*wi[0] );
  ws[1] = fabs( ww - deps*wi[1] );

  /* numerically approximate delta function, distribute between nearest array elements */
  if (wi[0] == wi[1]) {
    p_n[ wi[0] ] = 1.0/deps;
  } else {
    p_n[ wi[0] ] = (1.0 - ws[0]/deps)/deps;
    p_n[ wi[1] ] = (1.0 - ws[1]/deps)/deps;
  }

  /* calculate first passage time distribution */
  for (ii = 0; ii < N_dt - 1; ii++) {

    /* set decision threshold location at previous time step*/
    bu_n = bu_np1;
    dbudt_n = dbudt_np1;
    bl_n = bl_np1;
    dbldt_n = dbldt_np1;
    s_n = bu_n - bl_n;
    dsdt_n = dbudt_n - dbldt_n;

    /* set time step */
    if (ii < N_dt_ini) {
      dt_ = dt_ini;
    } else if ((ii >= N_dt_ini) && (ii < N_dt_ini + N_dt_scaleup)) {
      dt_ = dt_ini + (ii + 1.0 - N_dt_ini)*(dt_base - dt_ini)/(N_dt_scaleup);
    } else {
      dt_ = dt_base;
    }

    /* decrease time step if decision thresholds decrease too rapidly */
    jj = 0;
    kk = 0;
    while ((jj < 1) && (kk < N_decmax)) {
      t_in = tt + dt_;
      bu_np1 = upper_threshold(phi, t_in);
      bl_np1 = lower_threshold(phi, t_in);

      t_in = tt - dt_;
      if (t_in <= 0.0) {
        t_in = dt_ini/10.0;
        bu_nm1 = upper_threshold(phi, t_in);
        bl_nm1 = lower_threshold(phi, t_in);
      } else {
        bu_nm1 = upper_threshold(phi, t_in);
        bl_nm1 = lower_threshold(phi, t_in);
      }

      s_np1 = bu_np1 - bl_np1;
      s_nm1 = bu_nm1 - bl_nm1;

      ds_ratio = 0.0;
      ds_ratio_np1 = fabs(s_n - s_np1)/s_n;
      ds_ratio_nm1 = fabs(s_n - s_nm1)/s_n;

      if (ds_ratio_np1 < ds_ratio_nm1) {
        ds_ratio = ds_ratio_nm1;
      } else {
        ds_ratio = ds_ratio_np1;
      }

      if (ds_ratio > ds_ratio_cutoff) {
        dt_ = dt_*dt_mod_scale;
      } else {
        jj = 1;
      }

      kk += 1;

    }

    /* calculate new decision thresholds */
    dt_ = modify_dt(phi, tt)*dt_;
    tt += dt_;

    bu_np1 = upper_threshold(phi, tt);
    dbudt_np1 = ( upper_threshold(phi, tt+delt) - upper_threshold(phi, tt-delt) ) / (2.0*delt);

    if (bu_np1 <= threshold_cutoff) {
      bu_np1 = threshold_cutoff;
      dbudt_np1 = 0.0;
    }

    bl_np1 = lower_threshold(phi, tt);
    dbldt_np1 = ( lower_threshold(phi, tt+delt) - lower_threshold(phi, tt-delt) ) / (2.0*delt);

    if (bl_np1 >= -threshold_cutoff) {
      bl_np1 = -threshold_cutoff;
      dbldt_np1 = 0.0;
    }

    s_np1 = bu_np1 - bl_np1;
    dsdt_np1 = dbudt_np1 - dbldt_np1;

    /* calculate drift and diffusion rates */
    weight = ts_cdf(phi, tt);

    if (ii == 0) {
      v_n = drift(phi, tt-dt_, weight);
      sigma_n = diffusion(phi, tt-dt_, weight);
      sigma2_n = sigma_n*sigma_n;
    } else {
      sigma_n = sigma_np1;
      sigma2_n = sigma2_np1;
      v_n = v_np1;
    }
    v_np1 = drift(phi, tt, weight);
    sigma_np1 = diffusion(phi, tt, weight);
    sigma2_np1 = sigma_np1*sigma_np1;

    /* invert matrix using tridiagonal matrix algorithm */
    for (jj = 0; jj < 2; jj++) {
      p_np1[jj] = 0.0;
    }

    // for (jj = 0; jj < 2; jj++) {
    //
    //   p_np1[jj] = 0.0;
    //
    //   if (ii == 0) {
    //     x_n = s_n*eps[jj] + bl_n;
    //     sigma_n[jj] = diffusion(phi, x_n, tt-dt_);
    //     sigma2_n[jj] = sigma_n[jj]*sigma_n[jj];
    //   } else {
    //     x_n = x_np1;
    //     sigma_n[jj] = sigma_np1[jj];
    //     sigma2_n[jj] = sigma2_np1[jj];
    //   }
    //
    //   x_np1 = s_np1*eps[jj] + bl_np1;
    //   sigma_np1[jj] = diffusion(phi, x_np1, tt);
    //   sigma2_np1[jj] = sigma_np1[jj]*sigma_np1[jj];
    //
    // }

    for (jj = 1; jj < N_deps-1; jj++) {

      p_np1[jj+1] = 0.0;

      alpha_n = dt_/(4.0*s_n*deps);
      alpha_np1 = dt_/(4.0*s_np1*deps);

      beta_n = eps[jj]*dsdt_n + dbldt_n;
      beta_np1 = eps[jj]*dsdt_np1 + dbldt_np1;

      gamma_n = alpha_n/(s_n*deps);
      gamma_np1 = alpha_np1/(s_np1*deps);

      AA[jj] = alpha_np1*beta_np1 - alpha_np1*v_np1 - gamma_np1*sigma2_np1;
      BB[jj] = 1.0 + 2.0*gamma_np1*sigma2_np1;
      CC[jj] = -alpha_np1*beta_np1 + alpha_np1*v_np1 - gamma_np1*sigma2_np1;

      DD[jj] = -alpha_n*beta_n + alpha_n*v_n + gamma_n*sigma2_n;
      EE[jj] = 1.0 - 2.0*gamma_n*sigma2_n;
      FF[jj] = alpha_n*beta_n - alpha_n*v_n + gamma_n*sigma2_n;

      p_ncn[jj] = DD[jj]*p_n[jj-1] + EE[jj]*p_n[jj] + FF[jj]*p_n[jj+1];

    }

    for (jj = 2; jj < N_deps - 1; jj++) {
      WW = AA[jj]/BB[jj-1];
      BB[jj] = BB[jj] - WW*CC[jj-1];
      p_ncn[jj] = p_ncn[jj] - WW*p_ncn[jj-1];
    }

    p_np1[N_deps-2] = p_ncn[N_deps-2]/BB[N_deps-2];
    p_n[N_deps-2] = p_np1[N_deps-2];

    for (jj = N_deps - 3; jj > 0; jj--) {
      p_np1[jj] = (p_ncn[jj] - CC[jj]*p_np1[jj+1])/BB[jj];
      p_n[jj] = p_np1[jj];
    }

    /* calculate probability of first passage time */
    p_fpt[0][ii+1] = tt + tnd;
    p_fpt[1][ii+1] = fabs( 0.5 * sigma2_np1 * (4.0*p_np1[N_deps-2] - p_np1[N_deps-3])/(2.0*deps*s_np1*s_ini) );
    p_fpt[2][ii+1] = fabs( 0.5 * sigma2_np1 * (4.0*p_np1[1] - p_np1[2])/(2.0*deps*s_np1*s_ini) );

    /* if p_fpt less than p_fpt_min, set to p_fpt_min */
    if (p_fpt[1][ii+1] <= p_fpt_min) {
      p_fpt[1][ii+1] = p_fpt_min;
    }
    if (p_fpt[2][ii+1] <= p_fpt_min) {
      p_fpt[2][ii+1] = p_fpt_min;
    }

    /* determine integrated probability */
    int_prob += p_fpt[1][ii+1]*dt_ + p_fpt[2][ii+1]*dt_;

    /* stop solver if time above rt_max or probability below p_fpt_min */
    if (p_fpt[0][ii+1] > rt_max) {
      g_prob = 0.5*gg*contamination_probability(phi, p_fpt[0][ii+1]);
      p_fpt[1][ii+1] = g_prob + (1.0 - gg)*p_fpt[1][ii+1];
      p_fpt[2][ii+1] = g_prob + (1.0 - gg)*p_fpt[2][ii+1];
      cdf[0][ii+1] = cdf[0][ii] + 0.5*(p_fpt[0][ii+1] - p_fpt[0][ii])*(p_fpt[1][ii+1] + p_fpt[1][ii]);
      cdf[1][ii+1] = cdf[1][ii] + 0.5*(p_fpt[0][ii+1] - p_fpt[0][ii])*(p_fpt[2][ii+1] + p_fpt[2][ii]);
      N_cut = ii + 1;
      break;
    } else if ((int_prob > int_prob_min) && (p_fpt[1][ii+1] <= p_fpt_min) && (p_fpt[2][ii+1] <= p_fpt_min)) {
      g_prob = 0.5*gg*contamination_probability(phi, p_fpt[0][ii+1]);
      p_fpt[1][ii+1] = g_prob + (1.0 - gg)*p_fpt[1][ii+1];
      p_fpt[2][ii+1] = g_prob + (1.0 - gg)*p_fpt[2][ii+1];
      cdf[0][ii+1] = cdf[0][ii] + 0.5*(p_fpt[0][ii+1] - p_fpt[0][ii])*(p_fpt[1][ii+1] + p_fpt[1][ii]);
      cdf[1][ii+1] = cdf[1][ii] + 0.5*(p_fpt[0][ii+1] - p_fpt[0][ii])*(p_fpt[2][ii+1] + p_fpt[2][ii]);
      N_cut = ii + 1;
      break;
    }

    /* calculate cdf including contamination distribution */
    g_prob = 0.5*gg*contamination_probability(phi, p_fpt[0][ii+1]);
    p_fpt[1][ii+1] = g_prob + (1.0 - gg)*p_fpt[1][ii+1];
    p_fpt[2][ii+1] = g_prob + (1.0 - gg)*p_fpt[2][ii+1];
    cdf[0][ii+1] = cdf[0][ii] + 0.5*(p_fpt[0][ii+1] - p_fpt[0][ii])*(p_fpt[1][ii+1] + p_fpt[1][ii]);
    cdf[1][ii+1] = cdf[1][ii] + 0.5*(p_fpt[0][ii+1] - p_fpt[0][ii])*(p_fpt[2][ii+1] + p_fpt[2][ii]);

  }

  /* determine size of rt vectors */
  int N_rtu = rtu.size();
  int N_rtl = rtl.size();

  /* initialize logCDF variable */
  double logCDF = 0.0;

  /* interpolate to find logCDF of data */
  kk = 0;
  for (ii = 0; ii < N_rtu; ii++){
    if ( (rtu[ii] >= p_fpt[0][N_cut]) ) {
      RCDFupp[ii] = cdf[0][N_cut];
      RlogCDFupp[ii] = log(RCDFupp[ii]);
      logCDF += RlogCDFupp[ii];
    } else if ( (rtu[ii] <= p_fpt[0][0]) ) {
      RCDFupp[ii] = p_fpt_min;
      RlogCDFupp[ii] = log(RCDFupp[ii]);
      logCDF += RlogCDFupp[ii];
    } else {
      for (jj = kk; jj < N_cut; jj++){
        if ( (rtu[ii] >= p_fpt[0][jj]) && (rtu[ii] < p_fpt[0][jj+1]) ) {
          RCDFupp[ii] = cdf[0][jj] + (rtu[ii] - p_fpt[0][jj]) * (cdf[0][jj+1] - cdf[0][jj]) / (p_fpt[0][jj+1] - p_fpt[0][jj]);
          RlogCDFupp[ii] = log(RCDFupp[ii]);
          logCDF += RlogCDFupp[ii];
          kk = jj;
          break;
        }
      }
    }
  }

  kk = 0;
  for (ii = 0; ii < N_rtl; ii++){
    if ( (rtl[ii] >= p_fpt[0][N_cut]) ) {
      RCDFlow[ii] = cdf[1][N_cut];
      RlogCDFlow[ii] = log(RCDFlow[ii]);
      logCDF += RlogCDFlow[ii];
    } else if ( (rtl[ii] <= p_fpt[0][0]) ) {
      RCDFlow[ii] = p_fpt_min;
      RlogCDFlow[ii] = log(RCDFlow[ii]);
      logCDF += RlogCDFlow[ii];
    } else {
      for (jj = kk; jj < N_cut; jj++){
        if ( (rtl[ii] >= p_fpt[0][jj]) && (rtl[ii] < p_fpt[0][jj+1]) ) {
          RCDFlow[ii] = cdf[1][jj] + (rtl[ii] - p_fpt[0][jj]) * (cdf[1][jj+1] - cdf[1][jj]) / (p_fpt[0][jj+1] - p_fpt[0][jj]);
          RlogCDFlow[ii] = log(RCDFlow[ii]);
          logCDF += RlogCDFlow[ii];
          kk = jj;
          break;
        }
      }
    }
  }

  /* write logCDF to R object */
  RsumlogCDF[0] = logCDF;

  return 0;

}


/* method used to draw random samples */
int Model_TW::rand(double *Rrt, double *phi) const {

  /* initialize loop indicies */
  int ii = 0, jj = 0, kk = 0;

  /* initialize remaining variables */
  double tnd = non_decision(phi); /* relative start point */
  double ww = relative_start(phi); /* relative start point */
  double vv = 0.0; /* drift rate */
  double DD = diffusion(phi, 0.0, 0.0); /* diffusion rate */
  double bu = upper_threshold(phi, 0.0); /* upper decision threshold */
  double bl = lower_threshold(phi, 0.0); /* lower decision threshold */
  double zz = bl + ww*(bu - bl); /* calculate absolute start point */
  double sqrtdt = sqrt(dt_); /* calculate square root of time step */
  double tt = 0.0; /* time */
  double x_old = 0.0;
  double xx = 0.0; /* accumulated evidence */
  // double u1 = 0.0, u2 = 0.0,
  double uu = 0.0; /* random number value */
  double t_old = 0.0;
  double rt = 0.0; /* first passage time for threshold crossing */
  double weight = 0.0;

  /* seed rng */
  // srand(seed);
  GetRNGstate();

  /* simulate model */
  for (ii = 0; ii < N; ii++) {

    /* reset time and space for next simulation */
    tt = 0.0;
    xx = zz;
    jj = 0;

    while ((jj < 1) && (tt <= tsim_max)) {

      /* update drift rate, diffusion rate, and decision thresholds */
      weight = ts_cdf(phi, tt);
      vv = drift(phi, tt, weight);
      DD = diffusion(phi, tt, weight);

      /* update time step */
      t_old = tt;
      tt += dt_;

      /* update decision threshold */
      bu = upper_threshold(phi, tt);
      bl = lower_threshold(phi, tt);

      /* draw random number for simulating diffusion */
      // u1 = unif_L();
      // u2 = unif_L();
      // u1 = std::rand()/(RAND_MAX*1.0);
      // u2 = std::rand()/(RAND_MAX*1.0);
      // uu = sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
      // uu = norm_rand();
      uu = Rf_rnorm(0.0, sqrtdt);

      /* calculate accumulated evidence */
      xx += dt_*vv + DD*uu;

      /* update time step and decision threhsolds */
      tt += dt_;
      bu = upper_threshold(phi, tt);
      bl = lower_threshold(phi, tt);

      /* check if accumulated evidence has crossed a decision threshold */
      if (xx >= bu) {
        rt = t_old + (bu - x_old)/(xx - x_old)*(tt - t_old);
        Rrt[ii] = tnd + rt;
        jj = 1;
      } else if (xx <= bl) {
        rt = t_old + (bl - x_old)/(xx - x_old)*(tt - t_old);
        Rrt[ii] = -tnd - rt;
        jj = 1;
      }

    }
  }
  PutRNGstate();


  return 0;

}



/* method used to construct grid for PDF */
int Model_TW::grid_pdf(double *Rrt, double *Rpdf_u, double *Rpdf_l, double *phi) const {

  /* initialize loop indicies */
  int ii = 0, jj = 0, kk = 0;

  /* initialize arrays */
  double p_fpt[3][N_dt]; /* first passage time probability */
  double p_n[N_deps_max]; /* probability of evidence x at time step n */
  double p_ncn[N_deps_max]; /* Crank-Nicolson modified probability of evidence x at time step n */
  double p_np1[N_deps_max]; /* Crank-Nicolson modified probability of evidence x at time step n */
  double eps[N_deps_max]; /* probability of evidence x at time step n plus 1 */
  double AA[N_deps_max], BB[N_deps_max], CC[N_deps_max], DD[N_deps_max], EE[N_deps_max], FF[N_deps_max]; /* Crank-Nicolson matrix arrays */
  double ws[2]; /* start point array */
  int wi[2]; /* index of start point arrays */

  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < N_dt; jj++) {
      p_fpt[ii][jj] = 0.0;
    }
  }

  for (ii = 0; ii < N_deps_max; ii++) {
    p_n[ii] = 0.0;
    p_ncn[ii] = 0.0;
    p_np1[ii] = 0.0;
    eps[ii] = 0.0;
    AA[ii] = 0.0;
    BB[ii] = 0.0;
    CC[ii] = 0.0;
    DD[ii] = 0.0;
    EE[ii] = 0.0;
    FF[ii] = 0.0;
  }

  wi[0] = 0.0;
  wi[1] = 0.0;
  ws[0] = 0.0;
  ws[1] = 0.0;

  /* initilize constants */
  double v_n = 0.0, v_np1 = 0.0; /* drift rate at time steps n and n plus 1 */
  double sigma_n = 0.0, sigma_np1 = 0.0;  /* diffusion rate at time steps n and n plus 1 */
  double sigma2_n = 0.0, sigma2_np1 = 0.0; /* squared diffusion rate at time steps n and n plus 1 */
  double deps = 1.0/(N_deps - 1.0);
  double dt_base = 0.0, dt_ini = 0.0, dt_ = 0.0;
  double t_in = 0.0;
  double bu_nm1 = 0.0, bu_n = 0.0, bu_np1 = 0.0;
  double bl_nm1 = 0.0, bl_n = 0.0, bl_np1 = 0.0;
  double s_nm1 = 0.0, s_n = 0.0, s_np1 = 0.0;
  double dbudt_n = 0.0, dbudt_np1 = 0.0;
  double dbldt_n = 0.0, dbldt_np1 = 0.0;
  double dsdt_n = 0.0, dsdt_np1 = 0.0;
  double tnd = non_decision(phi);
  double ww = relative_start(phi);
  double tt = 0.0;
  double ds_ratio = 0.0, ds_ratio_np1 = 0.0, ds_ratio_nm1 = 0.0;
  double x_n = 0.0, x_np1 = 0.0;
  double alpha_n = 0.0, alpha_np1 = 0.0;
  double beta_n = 0.0, beta_np1 = 0.0;
  double gamma_n = 0.0, gamma_np1 = 0.0;
  double WW = 0.0;
  double s_ini = upper_threshold(phi, 0.0) - lower_threshold(phi, 0.0);
  double int_prob = 0.0;
  double gg = contamination_strength(phi);
  double g_prob = 0.0;
  double weight = 0.0;

  /* set value of first element for first passage time array, output to files */
  p_fpt[0][0] = tnd;
  g_prob = 0.5*gg*contamination_probability(phi, p_fpt[0][0]);
  p_fpt[1][0] = g_prob;
  p_fpt[2][0] = g_prob;
  Rrt[0] = p_fpt[0][0];
  Rpdf_u[0] = p_fpt[1][0];
  Rpdf_l[0] = p_fpt[2][0];

  /* create array of spatial locations */
  eps[0] = 0.0;
  for (ii = 0; ii < N_deps; ii++) {
    eps[ii] = ii*deps;
  }
  eps[N_deps-1] = 1.0;

  /* calculate time step sizes */
  dt_base = approx_dt(phi, dt_scale);
  if (dt_base > dt_max) {
    dt_base = dt_max;
  }
  dt_ini = dt_ini_scale*dt_base;
  dt_ = dt_ini;

  /* calculate initial threshold locations */
  t_in = dt_ini/10.0;
  bu_np1 = upper_threshold(phi, t_in);
  bl_np1 = lower_threshold(phi, t_in);

  /* calculate initial threshold derivatives */
  dbudt_np1 = ( upper_threshold(phi, t_in+delt) - upper_threshold(phi, t_in-delt) ) / (2.0*delt);
  dbldt_np1 = ( lower_threshold(phi, t_in+delt) - lower_threshold(phi, t_in-delt) ) / (2.0*delt);

  /* calculate initial threhsold separation and derivative */
  s_np1 = bu_np1 - bl_np1;
  dsdt_np1 = dbudt_np1 - dbldt_np1;

  /* calculate probability around start point */
  /* determine nearest array elements to delta distribution */
  wi[0] = int( ceil(ww/deps) );
  wi[1] = int( floor(ww/deps) );

  /* determine distance between delta distribution and nearest array elements */
  ws[0] = fabs( ww - deps*wi[0] );
  ws[1] = fabs( ww - deps*wi[1] );

  /* numerically approximate delta function, distribute between nearest array elements */
  if (wi[0] == wi[1]) {
    p_n[ wi[0] ] = 1.0/deps;
  } else {
    p_n[ wi[0] ] = (1.0 - ws[0]/deps)/deps;
    p_n[ wi[1] ] = (1.0 - ws[1]/deps)/deps;
  }

  /* calculate first passage time distribution */
  for (ii = 0; ii < N_dt - 1; ii++) {

    /* set decision threshold location at previous time step*/
    bu_n = bu_np1;
    dbudt_n = dbudt_np1;
    bl_n = bl_np1;
    dbldt_n = dbldt_np1;
    s_n = bu_n - bl_n;
    dsdt_n = dbudt_n - dbldt_n;

    /* set time step */
    if (ii < N_dt_ini) {
      dt_ = dt_ini;
    } else if ((ii >= N_dt_ini) && (ii < N_dt_ini + N_dt_scaleup)) {
      dt_ = dt_ini + (ii + 1.0 - N_dt_ini)*(dt_base - dt_ini)/(N_dt_scaleup);
    } else {
      dt_ = dt_base;
    }

    /* decrease time step if decision thresholds decrease too rapidly */
    jj = 0;
    kk = 0;
    while ((jj < 1) && (kk < N_decmax)) {
      t_in = tt + dt_;
      bu_np1 = upper_threshold(phi, t_in);
      bl_np1 = lower_threshold(phi, t_in);

      t_in = tt - dt_;
      if (t_in <= 0.0) {
        t_in = dt_ini/10.0;
        bu_nm1 = upper_threshold(phi, t_in);
        bl_nm1 = lower_threshold(phi, t_in);
      } else {
        bu_nm1 = upper_threshold(phi, t_in);
        bl_nm1 = lower_threshold(phi, t_in);
      }

      s_np1 = bu_np1 - bl_np1;
      s_nm1 = bu_nm1 - bl_nm1;

      ds_ratio = 0.0;
      ds_ratio_np1 = fabs(s_n - s_np1)/s_n;
      ds_ratio_nm1 = fabs(s_n - s_nm1)/s_n;

      if (ds_ratio_np1 < ds_ratio_nm1) {
        ds_ratio = ds_ratio_nm1;
      } else {
        ds_ratio = ds_ratio_np1;
      }

      if (ds_ratio > ds_ratio_cutoff) {
        dt_ = dt_*dt_mod_scale;
      } else {
        jj = 1;
      }

      kk += 1;

    }

    /* calculate new decision thresholds */
    dt_ = modify_dt(phi, tt)*dt_;
    tt += dt_;

    bu_np1 = upper_threshold(phi, tt);
    dbudt_np1 = ( upper_threshold(phi, tt+delt) - upper_threshold(phi, tt-delt) ) / (2.0*delt);

    if (bu_np1 <= threshold_cutoff) {
      bu_np1 = threshold_cutoff;
      dbudt_np1 = 0.0;
    }

    bl_np1 = lower_threshold(phi, tt);
    dbldt_np1 = ( lower_threshold(phi, tt+delt) - lower_threshold(phi, tt-delt) ) / (2.0*delt);

    if (bl_np1 >= -threshold_cutoff) {
      bl_np1 = -threshold_cutoff;
      dbldt_np1 = 0.0;
    }

    s_np1 = bu_np1 - bl_np1;
    dsdt_np1 = dbudt_np1 - dbldt_np1;

    /* calculate drift and diffusion rates */
    weight = ts_cdf(phi, tt);

    if (ii == 0) {
      v_n = drift(phi, tt-dt_, weight);
      sigma_n = diffusion(phi, tt-dt_, weight);
      sigma2_n = sigma_n*sigma_n;
    } else {
      sigma_n = sigma_np1;
      sigma2_n = sigma2_np1;
      v_n = v_np1;
    }
    v_np1 = drift(phi, tt, weight);
    sigma_np1 = diffusion(phi, tt, weight);
    sigma2_np1 = sigma_np1*sigma_np1;

    /* invert matrix using tridiagonal matrix algorithm */
    for (jj = 0; jj < 2; jj++) {
      p_np1[jj] = 0.0;
    }

    // for (jj = 0; jj < 2; jj++) {
    //
    //   p_np1[jj] = 0.0;
    //
    //   if (ii == 0) {
    //     x_n = s_n*eps[jj] + bl_n;
    //     sigma_n[jj] = diffusion(phi, x_n, tt-dt_);
    //     sigma2_n[jj] = sigma_n[jj]*sigma_n[jj];
    //   } else {
    //     x_n = x_np1;
    //     sigma_n[jj] = sigma_np1[jj];
    //     sigma2_n[jj] = sigma2_np1[jj];
    //   }
    //
    //   x_np1 = s_np1*eps[jj] + bl_np1;
    //   sigma_np1[jj] = diffusion(phi, x_np1, tt);
    //   sigma2_np1[jj] = sigma_np1[jj]*sigma_np1[jj];
    //
    // }

    for (jj = 1; jj < N_deps-1; jj++) {

      p_np1[jj+1] = 0.0;

      alpha_n = dt_/(4.0*s_n*deps);
      alpha_np1 = dt_/(4.0*s_np1*deps);

      beta_n = eps[jj]*dsdt_n + dbldt_n;
      beta_np1 = eps[jj]*dsdt_np1 + dbldt_np1;

      gamma_n = alpha_n/(s_n*deps);
      gamma_np1 = alpha_np1/(s_np1*deps);

      AA[jj] = alpha_np1*beta_np1 - alpha_np1*v_np1 - gamma_np1*sigma2_np1;
      BB[jj] = 1.0 + 2.0*gamma_np1*sigma2_np1;
      CC[jj] = -alpha_np1*beta_np1 + alpha_np1*v_np1 - gamma_np1*sigma2_np1;

      DD[jj] = -alpha_n*beta_n + alpha_n*v_n + gamma_n*sigma2_n;
      EE[jj] = 1.0 - 2.0*gamma_n*sigma2_n;
      FF[jj] = alpha_n*beta_n - alpha_n*v_n + gamma_n*sigma2_n;

      p_ncn[jj] = DD[jj]*p_n[jj-1] + EE[jj]*p_n[jj] + FF[jj]*p_n[jj+1];

    }

    for (jj = 2; jj < N_deps - 1; jj++) {
      WW = AA[jj]/BB[jj-1];
      BB[jj] = BB[jj] - WW*CC[jj-1];
      p_ncn[jj] = p_ncn[jj] - WW*p_ncn[jj-1];
    }

    p_np1[N_deps-2] = p_ncn[N_deps-2]/BB[N_deps-2];
    p_n[N_deps-2] = p_np1[N_deps-2];

    for (jj = N_deps - 3; jj > 0; jj--) {
      p_np1[jj] = (p_ncn[jj] - CC[jj]*p_np1[jj+1])/BB[jj];
      p_n[jj] = p_np1[jj];
    }

    /* calculate probability of first passage time */
    p_fpt[0][ii+1] = tt + tnd;
    p_fpt[1][ii+1] = fabs( 0.5 * sigma2_np1 * (4.0*p_np1[N_deps-2] - p_np1[N_deps-3])/(2.0*deps*s_np1*s_ini) );
    p_fpt[2][ii+1] = fabs( 0.5 * sigma2_np1 * (4.0*p_np1[1] - p_np1[2])/(2.0*deps*s_np1*s_ini) );

    /* if p_fpt less than p_fpt_min, set to p_fpt_min */
    if (p_fpt[1][ii+1] <= p_fpt_min) {
      p_fpt[1][ii+1] = p_fpt_min;
    }
    if (p_fpt[2][ii+1] <= p_fpt_min) {
      p_fpt[2][ii+1] = p_fpt_min;
    }

    /* determine integrated probability */
    int_prob += p_fpt[1][ii+1]*dt_ + p_fpt[2][ii+1]*dt_;

    /* stop solver if time above rt_max or probability below p_fpt_min */
    if (p_fpt[0][ii+1] > rt_max) {
      g_prob = 0.5*gg*contamination_probability(phi, p_fpt[0][ii+1]);
      p_fpt[1][ii+1] = g_prob + (1.0 - gg)*p_fpt[1][ii+1];
      p_fpt[2][ii+1] = g_prob + (1.0 - gg)*p_fpt[2][ii+1];
      Rrt[ii+1] = p_fpt[0][ii+1];
      Rpdf_u[ii+1] = p_fpt[1][ii+1];
      Rpdf_l[ii+1] = p_fpt[2][ii+1];
      break;
    } else if ((int_prob > int_prob_min) && (p_fpt[1][ii+1] <= p_fpt_min) && (p_fpt[2][ii+1] <= p_fpt_min)) {
      g_prob = 0.5*gg*contamination_probability(phi, p_fpt[0][ii+1]);
      p_fpt[1][ii+1] = g_prob + (1.0 - gg)*p_fpt[1][ii+1];
      p_fpt[2][ii+1] = g_prob + (1.0 - gg)*p_fpt[2][ii+1];
      Rrt[ii+1] = p_fpt[0][ii+1];
      Rpdf_u[ii+1] = p_fpt[1][ii+1];
      Rpdf_l[ii+1] = p_fpt[2][ii+1];
      break;
    }

    /* calculate pdf including contamination distribution */
    g_prob = 0.5*gg*contamination_probability(phi, p_fpt[0][ii+1]);
    p_fpt[1][ii+1] = g_prob + (1.0 - gg)*p_fpt[1][ii+1];
    p_fpt[2][ii+1] = g_prob + (1.0 - gg)*p_fpt[2][ii+1];

    Rrt[ii+1] = p_fpt[0][ii+1];
    Rpdf_u[ii+1] = p_fpt[1][ii+1];
    Rpdf_l[ii+1] = p_fpt[2][ii+1];

  }

  return 0;

}
