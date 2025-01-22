/* file: model_functions.cpp
 Functions for defining model.
 Author: Mathew Murrow and Raphael Hartmann
 Date: Jan 13, 2025 */

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

#include "Model_7P.h"


/* Model_7P Constructor */
Model_7P::Model_7P() {}

/* Model_7P Destructor */
Model_7P::~Model_7P() {}



/* method which sets time step for likelihood generating function */
double Model_7P::approx_dt(double* phi, double dt_scale) const {

  /* initialize variables */
  // int seed = 90210; /* set random number generator seed */
  int ii = 0, jj = 0; /* initialize loop indicies */
  double ww = relative_start_average(phi); /* relative start point */
  double vv = drift_average(phi); /* drift rate */
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
      // vv = drift(phi, xx, tt);
      // DD = diffusion(phi, xx, tt);
      // bu = upper_threshold(phi, tt);
      // bl = lower_threshold(phi, tt);

      /* draw random number for simulating diffusion */
      // uu = -1 +  2 * ( rand() % 2 );
      uu = unif_L() >= 0.5 ? 1.0 : -1.0;

      /* calculate accumulated evidence */
      xx += dt_sims*vv + sqrtdt*DD*uu;

      /* check if accumulated evidence has crossed a decision threshold */
      if ((xx >= bu) || (xx <= bl)) {
        fpt += tt;
        jj = 1;
      }

    }
  }
  PutRNGstate();

  /* if no threshold crossing, set time step to max */
  if (fpt == 0.0) {
    fpt = N_sims*t_max;
  }

  /* calculate estimated time step for likelihood generation function */
  dt_ = dt_scale*fpt/N_sims;

  /* output time step for likelihood generation function */
  return dt_;

}

/* ---------------------------------------------------- */
/* functions for calculating across-trial variabilities */
/* ---------------------------------------------------- */

/* pdf for standard normal distribution */
double Model_7P::snd_pdf(double x) const {
  double pdf = 1.0/sqrt(2.0*M_PI) * exp(-0.5*x*x);
  return pdf;
}

/* cdf for standard normal distribution */
double Model_7P::snd_cdf(double x) const {
  double cdf = 0.5*(1.0 + erf(x/sqrt(2)));
  return cdf;
}

/* function for interpolating an array */
double Model_7P::interp(double t, double t0, double t1, double y0, double y1) const {
  double yy = (y0*(t1 - t) +  y1*(t - t0)) / (t1 - t0);
  return yy;
}

/* uniform non-decision time distribution */
double Model_7P::tnd_uni(double t, double a, double b) const {
  double p_tnd = 0.0;
  if ((t >= a) && (t <= b)) {
    p_tnd = 1.0/(b - a);
  }
  return p_tnd;
}

/* truncated normal non-decision time distribution */
double Model_7P::tnd_tn(double* phi, double t) const {
  double mu = non_decision_average(phi);
  double sig = non_decision_width(phi);
  double p_tnd = (1.0/sig) * snd_pdf((t - mu)/sig) / ( 1.0 - snd_cdf(-mu/sig) );
  return p_tnd;
}

/* inverse gaussian (Wald) non-decision time distribution */
double Model_7P::tnd_ig(double* phi, double t) const {
  double mu = non_decision_average(phi);
  double lamb = non_decision_width(phi);
  double p_tnd = sqrt( lamb/(2.0*M_PI*t*t*t) ) * exp( -lamb*(t - mu)*(t - mu)/(2.0*mu*mu*t) );
  return p_tnd;
}

/* uniform relative start point distribution */
double Model_7P::w_uni(double eps, double a, double b) const {
  double p_w = 0.0;
  if ((eps >= a) && (eps <= b)) {
    p_w = 1.0/(b - a);
  }
  return p_w;
}

/* truncated normal relative start point distribution */
double Model_7P::w_tn(double* phi, double eps, double a, double b) const {
  double mu = relative_start_average(phi);
  double sig = relative_start_width(phi);
  double p_w = (1.0/sig) * snd_pdf((eps - mu)/sig) / ( snd_cdf((b - mu)/sig) - snd_cdf((a - mu)/sig) );
  return p_w;
}

/* function for uniform drift rate distribution */
double Model_7P::v_uni(double* phi, double v) const {
  double v_mu = drift_average(phi);
  double v_sig = drift_width(phi);
  double v_a = v_mu - v_sig;
  double v_b = v_mu + v_sig;
  double p_v = 1.0/(v_b - v_a);
  return p_v;
}

/* function for normal drift rate distribution */
double Model_7P::v_norm(double* phi, double v) const {
  double v_mu = drift_average(phi);
  double v_sig = drift_width(phi);
  double p_v = 1.0/sqrt(2.0*M_PI*v_sig*v_sig) * exp( -(v - v_mu)*(v - v_mu) / (2.0*v_sig*v_sig) );
  return p_v;
}


/* method for PDF */
int Model_7P::pdf(double *RsumlogPDF, double *RPDFlow, double *RPDFupp, double *RlogPDFlow, double *RlogPDFupp, std::vector<double> rtl, std::vector<double> rtu, double *phi) const {

  /* initialize loop indicies */
  int ii = 0, jj = 0, kk = 0, ll = 0, mm = 0, nn = 0;

  /* initialize arrays */
  double p_fpt[3][N_dt]; /* first passage time probability */
  double p_atv[2][N_dt]; /* first passage time probability for pdf with across trial variability in drift rate */
  double p_con[3][N_dt]; /* first passage time probability for pdf convolved with non-decision time distribution */
  double p_ini[N_deps_max]; /* initial probability distribution */
  double p_n[N_deps_max]; /* probability of evidence x at time step n */
  double p_ncn[N_deps_max]; /* Crank-Nicolson modified probability of evidence x at time step n */
  double p_np1[N_deps_max]; /* Crank-Nicolson modified probability of evidence x at time step n */
  double AA[N_deps_max], BB[N_deps_max], CC[N_deps_max], DD[N_deps_max], EE[N_deps_max], FF[N_deps_max]; /* Crank-Nicolson matrix arrays */
  double ws[2]; /* start point array */
  int wi[2]; /* index of start point arrays */

  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < N_dt; jj++) {
      p_fpt[ii][jj] = 0.0;
      p_con[ii][jj] = 0.0;
    }
  }

  for (ii = 0; ii < 2; ii++) {
    for (jj = 0; jj < N_dt; jj++) {
      p_atv[ii][jj] = 0.0;
    }
  }

  for (ii = 0; ii < N_deps_max; ii++) {
    p_ini[ii] = 0.0;
    p_n[ii] = 0.0;
    p_ncn[ii] = 0.0;
    p_np1[ii] = 0.0;
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
  double v_avg = drift_average(phi);
  double v_sd = drift_width(phi);
  double vv = 0.0;
  double dv = 0.0;
  double sigma = diffusion(phi, 0.0, 0.0);
  double sigma2 = sigma*sigma;
  double deps = 1.0/(N_deps - 1.0);
  double eps = 0.0;
  double p_ini_int = 0.0;
  double dt_base = 0.0, dt_ini = 0.0, dt_ = 0.0;
  double t_in = 0.0;
  double bu = upper_threshold(phi, 0.0), bl = lower_threshold(phi, 0.0);
  double ww = 0.0;
  double w_avg = relative_start_average(phi);
  double w_sig = relative_start_width(phi);
  double w_a = 0.0;
  double w_b = 0.0;
  double tt = 0.0;
  double alpha = 0.0;
  double gamma = 0.0;
  double WW = 0.0;
  double s_ini = bu - bl;
  double int_prob = 0.0;
  double gg = contamination_strength(phi);
  double g_prob = 0.0;
  double tnd_mu = non_decision_average(phi);
  double tnd_sig = non_decision_width(phi);
  double tnd_a = 0.0;
  double tnd_b = 0.0;
  double dcon = 0.0;
  double tau0 = 0.0;
  double tau1 = 0.0;
  double dtau = 0.0;
  double tmtau0 = 0.0;
  double tmtau1 = 0.0;
  double p_tnd = 0.0;
  int N_cut = 0;
  int N_print = 0;

  /* calculate drift rate step */
  if (v_dist == 0) {
    dv = ( (v_avg + v_sd) - (v_avg - v_sd) ) / (N_dv - 1.0);
  } else if (v_dist == 1) {
    dv = ( (v_avg + v_range*v_sd) - (v_avg - v_range*v_sd) ) / (N_dv - 1.0);
  }

  /* ------------------------------------------- */
  /* calculate relative start point distribution */
  /* ------------------------------------------- */

  /* delta relative start distribution */
  if (w_dist == 0) {

    ww = w_avg;

    /* calculate probability around start point, then determine nearest array elements to delta distribution */
    wi[0] = int( ceil(ww/deps) );
    wi[1] = int( floor(ww/deps) );

    /* determine distance between delta distribution and nearest array elements */
    ws[0] = fabs( ww - deps*wi[0] );
    ws[1] = fabs( ww - deps*wi[1] );

    /* numerically approximate delta function, distribute between nearest array elements */
    if (wi[0] == wi[1]) {
      p_ini[ wi[0] ] = 1.0/deps;
    } else {
      p_ini[ wi[0] ] = (1.0 - ws[0]/deps)/deps;
      p_ini[ wi[1] ] = (1.0 - ws[1]/deps)/deps;
    }

    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      p_ini_int += 0.5*deps*(p_ini[ii-1] + p_ini[ii]);
    }

    for (ii = 0; ii < N_deps; ii++) {
      p_ini[ii] = (1.0/p_ini_int)*p_ini[ii];
    }

  }

  /* uniform relative start distribution */
  if (w_dist == 1) {

    w_a = std::max(w_avg - w_sig, w_min);
    w_b = std::min(w_avg + w_sig, w_max);

    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      if ((w_a >= eps - deps) && (w_a < eps)) {
        if (fabs(w_a - (eps - deps)) <= fabs(eps - w_a)) {
          w_a = eps - deps;
        } else {
          w_a = eps;
        }
        break;
      }
    }

    for (ii = N_deps; ii > 0; ii--) {
      eps = ii*deps;
      if ((w_b >= eps - deps) && (w_b < eps)) {
        if (fabs(w_b - (eps - deps)) <= fabs(eps - w_b)) {
          w_b = eps - deps;
        } else {
          w_b = eps;
        }
        break;
      }
    }

    p_ini[0] = w_uni(eps, w_a, w_b);
    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      p_ini[ii] = w_uni(eps, w_a, w_b);
      p_ini_int += 0.5*deps*(p_ini[ii-1] + p_ini[ii]);
    }

    for (ii = 0; ii < N_deps; ii++) {
      p_ini[ii] = 1.0/p_ini_int*p_ini[ii];
    }

  }

  /* truncated normal relative start distribution */
  if (w_dist == 2) {

    w_a = w_min;
    w_b = w_max;

    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      if ((w_a >= eps - deps) && (w_a < eps)) {
        if (fabs(w_a - (eps - deps)) <= fabs(eps - w_a)) {
          w_a = eps - deps;
        } else {
          w_a = eps;
        }
        break;
      }
    }

    for (ii = N_deps; ii > 0; ii--) {
      eps = ii*deps;
      if ((w_b >= eps - deps) && (w_b < eps)) {
        if (fabs(w_b - (eps - deps)) <= fabs(eps - w_b)) {
          w_b = eps - deps;
        } else {
          w_b = eps;
        }
        break;
      }
    }

    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      if ((eps >= w_a) && (eps <= w_b)) {
        p_ini[ii] = w_tn(phi, eps, w_min, w_max);
      }
      p_ini_int += 0.5*deps*(p_ini[ii-1] + p_ini[ii]);
    }

    for (ii = 0; ii < N_deps; ii++) {
      p_ini[ii] = (1.0/p_ini_int)*p_ini[ii];
    }

  }

  /* -------------------------------------------------------------------------------------------- */
  /* calculate first passage time distribution N_dv times for drift rate across trial variability */
  /* -------------------------------------------------------------------------------------------- */

  /* calculate time step sizes */
  dt_base = approx_dt(phi, dt_scale);
  if (dt_base > dt_max) {
    dt_base = dt_max;
  }
  dt_ini = dt_ini_scale*dt_base;
  dt_ = dt_ini;

  /* set non-decision width to dt_base if smaller */
  tnd_sig = std::max(tnd_sig, dt_base);
  tnd_a = tnd_mu - tnd_sig;
  tnd_b = tnd_mu + tnd_sig;

  for (nn = 0; nn < N_dv; nn++) {

    /* reset arrays */
    for (ii = 1; ii < 3; ii++) {
      for (jj = 0; jj < N_dt; jj++) {
        p_fpt[ii][jj] = 0.0;
      }
    }

    for (ii = 0; ii < N_deps_max; ii++) {
      p_n[ii] = p_ini[ii];
      p_ncn[ii] = 0.0;
      p_np1[ii] = 0.0;
      AA[ii] = 0.0;
      BB[ii] = 0.0;
      CC[ii] = 0.0;
      DD[ii] = 0.0;
      EE[ii] = 0.0;
      FF[ii] = 0.0;
    }

    /* reset constants */
    if (v_dist == 0){
      vv = v_avg - v_sd + nn*dv;
    } else if (v_dist == 1) {
      vv = v_avg - v_range*v_sd + nn*dv;
    }

    tt = 0.0;
    alpha = 0.0;
    gamma = 0.0;
    WW = 0.0;
    int_prob = 0.0;
    g_prob = 0.0;

    /* calculate first passage time distribution */
    for (ii = 0; ii < N_dt - 1; ii++) {

      /* set time step */
      if (ii < N_dt_ini) {
        dt_ = dt_ini;
      } else if ((ii >= N_dt_ini) && (ii < N_dt_ini + N_dt_scaleup)) {
        dt_ = dt_ini + (ii + 1.0 - N_dt_ini)*(dt_base - dt_ini)/(N_dt_scaleup);
      } else {
        dt_ = dt_base;
      }
      dt_ = modify_dt(phi, tt)*dt_;
      tt += dt_;

      /* invert matrix using tridiagonal matrix algorithm */
      p_np1[0] = 0.0;
      p_np1[1] = 0.0;

      for (jj = 1; jj < N_deps-1; jj++) {

        p_np1[jj+1] = 0.0;

        alpha = dt_/(4.0*s_ini*deps);
        gamma = alpha/(s_ini*deps);

        AA[jj] = -alpha*vv - gamma*sigma2;
        BB[jj] = 1.0 + 2.0*gamma*sigma2;
        CC[jj] = alpha*vv - gamma*sigma2;

        DD[jj] = alpha*vv + gamma*sigma2;
        EE[jj] = 1.0 - 2.0*gamma*sigma2;
        FF[jj] = -alpha*vv + gamma*sigma2;

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
      p_fpt[0][ii+1] = tt;
      p_fpt[1][ii+1] = fabs( 0.5 * sigma2 * (4.0*p_np1[N_deps-2] - p_np1[N_deps-3])/(2.0*deps*s_ini*s_ini) );
      p_fpt[2][ii+1] = fabs( 0.5 * sigma2 * (4.0*p_np1[1] - p_np1[2])/(2.0*deps*s_ini*s_ini) );

      /* determine integrated probability */
      int_prob += p_fpt[1][ii+1]*dt_ + p_fpt[2][ii+1]*dt_;

      if (v_dist == 0) {
        if ( (nn == 0) || (nn == N_dv - 1) ) {
          p_atv[0][ii+1] += dv*v_uni(phi, vv)*p_fpt[1][ii+1]/2.0;
          p_atv[1][ii+1] += dv*v_uni(phi, vv)*p_fpt[2][ii+1]/2.0;
        } else {
          p_atv[0][ii+1] += dv*v_uni(phi, vv)*p_fpt[1][ii+1];
          p_atv[1][ii+1] += dv*v_uni(phi, vv)*p_fpt[2][ii+1];
        }
      } else if (v_dist == 1) {
        if ( (nn == 0) || (nn == N_dv - 1) ) {
          p_atv[0][ii+1] += dv*v_norm(phi, vv)*p_fpt[1][ii+1]/2.0;
          p_atv[1][ii+1] += dv*v_norm(phi, vv)*p_fpt[2][ii+1]/2.0;
        } else {
          p_atv[0][ii+1] += dv*v_norm(phi, vv)*p_fpt[1][ii+1];
          p_atv[1][ii+1] += dv*v_norm(phi, vv)*p_fpt[2][ii+1];
        }
      }

      /* determine index of highest time point */
      N_cut = std::max(ii+1, N_cut);

      /* stop solver if time above rt_max or probability below p_fpt_min */
      if (p_fpt[0][ii+1] > rt_max) {
        break;
      } else if ((int_prob > int_prob_min) && (p_fpt[1][ii+1] <= p_fpt_min) && (p_fpt[2][ii+1] <= p_fpt_min)) {
        break;
      }

    }

  }

  /* ------------------------------------------ */
  /* convolve pdf and non-decision distribution */
  /* ------------------------------------------ */

  /* delta non-decision time distribution */
  if (tnd_dist == 0) {

    dcon = tnd_mu/(N_con - 1.0);
    for (ii = 0; ii <= N_cut + N_con - 1; ii++) {
      if (ii < N_con) {
        p_con[0][ii] = ii*dcon;
      } else {
        p_con[0][ii] = p_fpt[0][ii-N_con+1] + tnd_mu;
      }
    }

    ll = 1;

    for (ii = 1; ii <= N_cut+N_con-1; ii++) {

      tmtau0 = p_con[0][ii] - tnd_mu;

      if (tmtau0 > 0.0) {
        for (kk = ll; kk <= N_cut; kk++) {
          if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
            p_con[1][ii] = interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
            p_con[2][ii] = interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
            ll = kk;
            break;
          }
        }
      }

    }
  }

  /* uniform non-decision time distribution */
  if (tnd_dist == 1) {

    tnd_b = std::min(tnd_b, p_fpt[0][N_cut]);
    dtau = (tnd_b - tnd_a)/(N_dtau - 1.0);

    dcon = tnd_a/(N_con - 1.0);
    for (ii = 0; ii <= N_cut + N_con - 1; ii++) {
      if (ii < N_con) {
        p_con[0][ii] = ii*dcon;
      } else {
        p_con[0][ii] = p_fpt[0][ii-N_con+1] + tnd_a;
      }
    }

    for (ii = 1; ii <= N_cut+N_con-1; ii++) {

      ll = N_cut;
      mm = N_cut;

      for (jj = 1; jj < N_dtau; jj++){

        tau0 = tnd_a + (jj-1.0)*dtau;
        tau1 = tnd_a + jj*dtau;
        tmtau0 = p_con[0][ii] - tau0;
        tmtau1 = p_con[0][ii] - tau1;

        if (tmtau0 > 0.0) {
          for (kk = ll; kk > 0; kk--) {
            if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
              p_tnd = tnd_uni(tau0, tnd_a, tnd_b);
              p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
              p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
              ll = kk;
              break;
            }
          }
        }

        if (tmtau1 > 0.0) {
          for (kk = mm; kk > 0; kk--) {
            if ((tmtau1 >= p_fpt[0][kk-1]) && (tmtau1 < p_fpt[0][kk])) {
              p_tnd = tnd_uni(tau1, tnd_a, tnd_b);
              p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
              p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
              mm = kk;
              break;
            }
          }
        }

      }
    }
  }

  /* truncated normal non-decision time distribution */
  if (tnd_dist == 2) {

    tnd_a = tnd_mu - tnd_range*tnd_sig;
    tnd_a = std::max(tnd_a, 0.0);
    tnd_b = tnd_mu + tnd_range*tnd_sig;
    tnd_b = std::min(tnd_b, p_fpt[0][N_cut]);
    dtau = (tnd_b - tnd_a)/(N_dtau - 1.0);

    if (tnd_a > 0.0) {

      dcon = (tnd_a)/(N_con - 1.0);
      for (ii = 0; ii <= N_cut + N_con - 1; ii++) {
        if (ii < N_con) {
          p_con[0][ii] = ii*dcon;
        } else {
          p_con[0][ii] = p_fpt[0][ii-N_con+1] + tnd_a;
        }
      }

      for (ii = 1; ii <= N_cut+N_con-1; ii++) {

        ll = N_cut;
        mm = N_cut;

        for (jj = 1; jj < N_dtau; jj++){

          tau0 = tnd_a + (jj-1.0)*dtau;
          tau1 = tnd_a + jj*dtau;
          tmtau0 = p_con[0][ii] - tau0;
          tmtau1 = p_con[0][ii] - tau1;

          if (tmtau0 > 0.0) {
            for (kk = ll; kk > 0; kk--) {
              if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
                p_tnd = tnd_tn(phi, tau0);
                p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
                p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
                ll = kk;
                break;
              }
            }
          }

          if (tmtau1 > 0.0) {
            for (kk = mm; kk > 0; kk--) {
              if ((tmtau1 >= p_fpt[0][kk-1]) && (tmtau1 < p_fpt[0][kk])) {
                p_tnd = tnd_tn(phi, tau1);
                p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
                p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
                mm = kk;
                break;
              }
            }
          }

        }
      }
    } else {

      for (ii = 0; ii <= N_cut; ii++) {
        p_con[0][ii] = p_fpt[0][ii];
      }

      for (ii = 1; ii <= N_cut; ii++) {

        ll = N_cut;
        mm = N_cut;

        for (jj = 1; jj < N_dtau; jj++){

          tau0 = (jj-1.0)*dtau;
          tau1 = jj*dtau;
          tmtau0 = p_con[0][ii] - tau0;
          tmtau1 = p_con[0][ii] - tau1;

          if (tmtau0 > 0.0) {
            for (kk = ll; kk > 0; kk--) {
              if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
                p_tnd = tnd_tn(phi, tau0);
                p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
                p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
                ll = kk;
                break;
              }
            }
          }

          if (tmtau1 > 0.0) {
            for (kk = mm; kk > 0; kk--) {
              if ((tmtau1 >= p_fpt[0][kk-1]) && (tmtau1 < p_fpt[0][kk])) {
                p_tnd = tnd_tn(phi, tau1);
                p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
                p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
                mm = kk;
                break;
              }
            }
          }

        }
      }
    }
  }

  /* inverse gaussian non-decision time distribution */
  if (tnd_dist == 3) {

    for (ii = 0; ii <= N_cut; ii++) {
      p_con[0][ii] = p_fpt[0][ii];
    }

    for (ii = 1; ii <= N_cut; ii++) {

      ll = N_cut;
      mm = N_cut;

      for (jj = 1; jj <= N_cut; jj++){

        tau0 = p_con[0][jj-1];
        tau1 = p_con[0][jj];
        dtau = tau1 - tau0;

        tmtau0 = p_con[0][ii] - tau0;
        tmtau1 = p_con[0][ii] - tau1;

        if (tmtau0 > 0.0) {
          for (kk = ll; kk > 0; kk--) {
            if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
              if (tau0 <= 0.0) {
                p_tnd = 0.0;
              } else {
                p_tnd = tnd_ig(phi, tau0);
              }
              p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
              p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
              ll = kk;
              break;
            }
          }
        }

        if (tmtau1 > 0.0) {
          for (kk = mm; kk > 0; kk--) {
            if ((tmtau1 >= p_fpt[0][kk-1]) && (tmtau1 < p_fpt[0][kk])) {
              if (tau1 <= 0.0) {
                p_tnd = 0.0;
              } else {
                p_tnd = tnd_ig(phi, tau1);
              }
              p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
              p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
              mm = kk;
              break;
            }
          }
        }

      }
    }
  }

  /* ----------------------- */
  /* calculate probabilities */
  /* ----------------------- */

  /* include "dummy" time steps into loop count */
  N_print = N_cut;
  if ((tnd_dist == 0) || (tnd_dist == 1)) {
    N_print = N_cut + N_con - 1;
  }
  if ((tnd_dist == 2) && (tnd_a > 0.0)) {
    N_print = N_cut + N_con - 1;
  }

  for (ii = 0; ii <= N_print; ii++) {

    /* include contamination distribution */
    g_prob = 0.5*gg*contamination_probability(phi, p_con[0][ii]);
    p_con[1][ii] = (1.0 - gg)*p_con[1][ii] + g_prob;
    p_con[2][ii] = (1.0 - gg)*p_con[2][ii] + g_prob;

    /* set small values to p_fpt_min */
    if (p_con[1][ii] < p_fpt_min) {
      p_con[1][ii] = p_fpt_min;
    }
    if (p_con[2][ii] < p_fpt_min) {
      p_con[2][ii] = p_fpt_min;
    }

  }

  /* --------------------------- */
  /* calculate logpdf of rt data */
  /* --------------------------- */

  /* determine size of rt vectors */
  int N_rtu = rtu.size();
  int N_rtl = rtl.size();

  /* initialize loglikelihood variable */
  double logPDF = 0.0;

  /* interpolate to find loglikelihood of data */
  kk = 0;
  for (ii = 0; ii < N_rtu; ii++){
    if ( (rtu[ii] > p_con[0][N_print]) || (rtu[ii] < p_con[0][0]) ) {
      RPDFupp[ii] = p_fpt_min;
      RlogPDFupp[ii] = log(RPDFupp[ii]);
      logPDF += RlogPDFupp[ii];
    } else {
      for (jj = kk; jj < N_print; jj++){
        if ( (rtu[ii] >= p_con[0][jj]) && (rtu[ii] < p_con[0][jj+1]) ) {
          RPDFupp[ii] = p_con[1][jj] + (rtu[ii] - p_con[0][jj]) * (p_con[1][jj+1] - p_con[1][jj]) / (p_con[0][jj+1] - p_con[0][jj]);
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
    if ( (rtl[ii] > p_con[0][N_print]) || (rtl[ii] < p_con[0][0]) ) {
      RPDFlow[ii] = p_fpt_min;
      RlogPDFlow[ii] = log(RPDFlow[ii]);
      logPDF += RlogPDFlow[ii];
    } else {
      for (jj = kk; jj < N_print; jj++){
        if ( (rtl[ii] >= p_con[0][jj]) && (rtl[ii] < p_con[0][jj+1]) ) {
          RPDFlow[ii] = p_con[2][jj] + (rtl[ii] - p_con[0][jj]) * (p_con[2][jj+1] - p_con[2][jj]) / (p_con[0][jj+1] - p_con[0][jj]);
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
int Model_7P::cdf(double *RsumlogCDF, double *RCDFlow, double *RCDFupp, double *RlogCDFlow, double *RlogCDFupp, std::vector<double> rtl, std::vector<double> rtu, double *phi) const {

  /* initialize loop indicies */
  int ii = 0, jj = 0, kk = 0, ll = 0, mm = 0, nn = 0;

  /* initialize arrays */
  double p_fpt[3][N_dt]; /* first passage time probability */
  double cdf[2][N_dt]; /* cumulative distribution function of first passage time probability */
  double p_atv[2][N_dt]; /* first passage time probability for pdf with across trial variability in drift rate */
  double p_con[3][N_dt]; /* first passage time probability for pdf convolved with non-decision time distribution */
  double p_ini[N_deps_max]; /* initial probability distribution */
  double p_n[N_deps_max]; /* probability of evidence x at time step n */
  double p_ncn[N_deps_max]; /* Crank-Nicolson modified probability of evidence x at time step n */
  double p_np1[N_deps_max]; /* Crank-Nicolson modified probability of evidence x at time step n */
  double AA[N_deps_max], BB[N_deps_max], CC[N_deps_max], DD[N_deps_max], EE[N_deps_max], FF[N_deps_max]; /* Crank-Nicolson matrix arrays */
  double ws[2]; /* start point array */
  int wi[2]; /* index of start point arrays */

  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < N_dt; jj++) {
      p_fpt[ii][jj] = 0.0;
      p_con[ii][jj] = 0.0;
    }
  }

  for (ii = 0; ii < 2; ii++) {
    for (jj = 0; jj < N_dt; jj++) {
      p_atv[ii][jj] = 0.0;
      cdf[ii][jj] = 0.0;
    }
  }

  for (ii = 0; ii < N_deps_max; ii++) {
    p_ini[ii] = 0.0;
    p_n[ii] = 0.0;
    p_ncn[ii] = 0.0;
    p_np1[ii] = 0.0;
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
  double v_avg = drift_average(phi);
  double v_sd = drift_width(phi);
  double vv = 0.0;
  double dv = 0.0;
  double sigma = diffusion(phi, 0.0, 0.0);
  double sigma2 = sigma*sigma;
  double deps = 1.0/(N_deps - 1.0);
  double eps = 0.0;
  double p_ini_int = 0.0;
  double dt_base = 0.0, dt_ini = 0.0, dt_ = 0.0;
  double t_in = 0.0;
  double bu = upper_threshold(phi, 0.0), bl = lower_threshold(phi, 0.0);
  double ww = 0.0;
  double w_avg = relative_start_average(phi);
  double w_sig = relative_start_width(phi);
  double w_a = 0.0;
  double w_b = 0.0;
  double tt = 0.0;
  double alpha = 0.0;
  double gamma = 0.0;
  double WW = 0.0;
  double s_ini = bu - bl;
  double int_prob = 0.0;
  double gg = contamination_strength(phi);
  double g_prob = 0.0;
  double tnd_mu = non_decision_average(phi);
  double tnd_sig = non_decision_width(phi);
  double tnd_a = tnd_mu - tnd_sig;
  double tnd_b = tnd_mu + tnd_sig;
  double dcon = 0.0;
  double tau0 = 0.0;
  double tau1 = 0.0;
  double dtau = 0.0;
  double tmtau0 = 0.0;
  double tmtau1 = 0.0;
  double p_tnd = 0.0;
  int N_cut = 0;
  int N_print = 0;

  /* calculate drift rate step */
  if (v_dist == 0) {
    dv = ( (v_avg + v_sd) - (v_avg - v_sd) ) / (N_dv - 1.0);
  } else if (v_dist == 1) {
    dv = ( (v_avg + v_range*v_sd) - (v_avg - v_range*v_sd) ) / (N_dv - 1.0);
  }

  /* ------------------------------------------- */
  /* calculate relative start point distribution */
  /* ------------------------------------------- */

  /* delta relative start distribution */
  if (w_dist == 0) {

    ww = w_avg;

    /* calculate probability around start point, then determine nearest array elements to delta distribution */
    wi[0] = int( ceil(ww/deps) );
    wi[1] = int( floor(ww/deps) );

    /* determine distance between delta distribution and nearest array elements */
    ws[0] = fabs( ww - deps*wi[0] );
    ws[1] = fabs( ww - deps*wi[1] );

    /* numerically approximate delta function, distribute between nearest array elements */
    if (wi[0] == wi[1]) {
      p_ini[ wi[0] ] = 1.0/deps;
    } else {
      p_ini[ wi[0] ] = (1.0 - ws[0]/deps)/deps;
      p_ini[ wi[1] ] = (1.0 - ws[1]/deps)/deps;
    }

    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      p_ini_int += 0.5*deps*(p_ini[ii-1] + p_ini[ii]);
    }

    for (ii = 0; ii < N_deps; ii++) {
      p_ini[ii] = (1.0/p_ini_int)*p_ini[ii];
    }

  }

  /* uniform relative start distribution */
  if (w_dist == 1) {

    w_a = std::max(w_avg - w_sig, w_min);
    w_b = std::min(w_avg + w_sig, w_max);

    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      if ((w_a >= eps - deps) && (w_a < eps)) {
        if (fabs(w_a - (eps - deps)) <= fabs(eps - w_a)) {
          w_a = eps - deps;
        } else {
          w_a = eps;
        }
        break;
      }
    }

    for (ii = N_deps; ii > 0; ii--) {
      eps = ii*deps;
      if ((w_b >= eps - deps) && (w_b < eps)) {
        if (fabs(w_b - (eps - deps)) <= fabs(eps - w_b)) {
          w_b = eps - deps;
        } else {
          w_b = eps;
        }
        break;
      }
    }

    p_ini[0] = w_uni(eps, w_a, w_b);
    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      p_ini[ii] = w_uni(eps, w_a, w_b);
      p_ini_int += 0.5*deps*(p_ini[ii-1] + p_ini[ii]);
    }

    for (ii = 0; ii < N_deps; ii++) {
      p_ini[ii] = 1.0/p_ini_int*p_ini[ii];
    }

  }

  /* truncated normal relative start distribution */
  if (w_dist == 2) {

    w_a = w_min;
    w_b = w_max;

    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      if ((w_a >= eps - deps) && (w_a < eps)) {
        if (fabs(w_a - (eps - deps)) <= fabs(eps - w_a)) {
          w_a = eps - deps;
        } else {
          w_a = eps;
        }
        break;
      }
    }

    for (ii = N_deps; ii > 0; ii--) {
      eps = ii*deps;
      if ((w_b >= eps - deps) && (w_b < eps)) {
        if (fabs(w_b - (eps - deps)) <= fabs(eps - w_b)) {
          w_b = eps - deps;
        } else {
          w_b = eps;
        }
        break;
      }
    }

    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      if ((eps >= w_a) && (eps <= w_b)) {
        p_ini[ii] = w_tn(phi, eps, w_min, w_max);
      }
      p_ini_int += 0.5*deps*(p_ini[ii-1] + p_ini[ii]);
    }

    for (ii = 0; ii < N_deps; ii++) {
      p_ini[ii] = (1.0/p_ini_int)*p_ini[ii];
    }

  }

  /* -------------------------------------------------------------------------------------------- */
  /* calculate first passage time distribution N_dv times for drift rate across trial variability */
  /* -------------------------------------------------------------------------------------------- */

  /* calculate time step sizes */
  dt_base = approx_dt(phi, dt_scale);
  if (dt_base > dt_max) {
    dt_base = dt_max;
  }
  dt_ini = dt_ini_scale*dt_base;
  dt_ = dt_ini;

  /* set non-decision width to dt_base if smaller */
  tnd_sig = std::max(tnd_sig, dt_base);

  for (nn = 0; nn < N_dv; nn++) {

    /* reset arrays */
    for (ii = 1; ii < 3; ii++) {
      for (jj = 0; jj < N_dt; jj++) {
        p_fpt[ii][jj] = 0.0;
      }
    }

    for (ii = 0; ii < N_deps_max; ii++) {
      p_n[ii] = p_ini[ii];
      p_ncn[ii] = 0.0;
      p_np1[ii] = 0.0;
      AA[ii] = 0.0;
      BB[ii] = 0.0;
      CC[ii] = 0.0;
      DD[ii] = 0.0;
      EE[ii] = 0.0;
      FF[ii] = 0.0;
    }

    /* reset constants */
    if (v_dist == 0){
      vv = v_avg - v_sd + nn*dv;
    } else if (v_dist == 1) {
      vv = v_avg - v_range*v_sd + nn*dv;
    }

    tt = 0.0;
    alpha = 0.0;
    gamma = 0.0;
    WW = 0.0;
    int_prob = 0.0;
    g_prob = 0.0;

    /* calculate first passage time distribution */
    for (ii = 0; ii < N_dt - 1; ii++) {

      /* set time step */
      if (ii < N_dt_ini) {
        dt_ = dt_ini;
      } else if ((ii >= N_dt_ini) && (ii < N_dt_ini + N_dt_scaleup)) {
        dt_ = dt_ini + (ii + 1.0 - N_dt_ini)*(dt_base - dt_ini)/(N_dt_scaleup);
      } else {
        dt_ = dt_base;
      }
      dt_ = modify_dt(phi, tt)*dt_;
      tt += dt_;

      /* invert matrix using tridiagonal matrix algorithm */
      p_np1[0] = 0.0;
      p_np1[1] = 0.0;

      for (jj = 1; jj < N_deps-1; jj++) {

        p_np1[jj+1] = 0.0;

        alpha = dt_/(4.0*s_ini*deps);
        gamma = alpha/(s_ini*deps);

        AA[jj] = -alpha*vv - gamma*sigma2;
        BB[jj] = 1.0 + 2.0*gamma*sigma2;
        CC[jj] = alpha*vv - gamma*sigma2;

        DD[jj] = alpha*vv + gamma*sigma2;
        EE[jj] = 1.0 - 2.0*gamma*sigma2;
        FF[jj] = -alpha*vv + gamma*sigma2;

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
      p_fpt[0][ii+1] = tt;
      p_fpt[1][ii+1] = fabs( 0.5 * sigma2 * (4.0*p_np1[N_deps-2] - p_np1[N_deps-3])/(2.0*deps*s_ini*s_ini) );
      p_fpt[2][ii+1] = fabs( 0.5 * sigma2 * (4.0*p_np1[1] - p_np1[2])/(2.0*deps*s_ini*s_ini) );

      /* determine integrated probability */
      int_prob += p_fpt[1][ii+1]*dt_ + p_fpt[2][ii+1]*dt_;

      if (v_dist == 0) {
        if ( (nn == 0) || (nn == N_dv - 1) ) {
          p_atv[0][ii+1] += dv*v_uni(phi, vv)*p_fpt[1][ii+1]/2.0;
          p_atv[1][ii+1] += dv*v_uni(phi, vv)*p_fpt[2][ii+1]/2.0;
        } else {
          p_atv[0][ii+1] += dv*v_uni(phi, vv)*p_fpt[1][ii+1];
          p_atv[1][ii+1] += dv*v_uni(phi, vv)*p_fpt[2][ii+1];
        }
      } else if (v_dist == 1) {
        if ( (nn == 0) || (nn == N_dv - 1) ) {
          p_atv[0][ii+1] += dv*v_norm(phi, vv)*p_fpt[1][ii+1]/2.0;
          p_atv[1][ii+1] += dv*v_norm(phi, vv)*p_fpt[2][ii+1]/2.0;
        } else {
          p_atv[0][ii+1] += dv*v_norm(phi, vv)*p_fpt[1][ii+1];
          p_atv[1][ii+1] += dv*v_norm(phi, vv)*p_fpt[2][ii+1];
        }
      }

      /* determine index of highest time point */
      N_cut = std::max(ii+1, N_cut);

      /* stop solver if time above rt_max or probability below p_fpt_min */
      if (p_fpt[0][ii+1] > rt_max) {
        break;
      } else if ((int_prob > int_prob_min) && (p_fpt[1][ii+1] <= p_fpt_min) && (p_fpt[2][ii+1] <= p_fpt_min)) {
        break;
      }

    }

  }

  /* ------------------------------------------ */
  /* convolve pdf and non-decision distribution */
  /* ------------------------------------------ */

  /* delta non-decision time distribution */
  if (tnd_dist == 0) {

    dcon = tnd_mu/(N_con - 1.0);
    for (ii = 0; ii <= N_cut + N_con - 1; ii++) {
      if (ii < N_con) {
        p_con[0][ii] = ii*dcon;
      } else {
        p_con[0][ii] = p_fpt[0][ii-N_con+1] + tnd_mu;
      }
    }

    ll = 1;

    for (ii = 1; ii <= N_cut+N_con-1; ii++) {

      tmtau0 = p_con[0][ii] - tnd_mu;

      if (tmtau0 > 0.0) {
        for (kk = ll; kk <= N_cut; kk++) {
          if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
            p_con[1][ii] = interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
            p_con[2][ii] = interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
            ll = kk;
            break;
          }
        }
      }

    }
  }

  /* uniform non-decision time distribution */
  if (tnd_dist == 1) {

    tnd_b = std::min(tnd_b, p_fpt[0][N_cut]);
    dtau = (tnd_b - tnd_a)/(N_dtau - 1.0);

    dcon = tnd_a/(N_con - 1.0);
    for (ii = 0; ii <= N_cut + N_con - 1; ii++) {
      if (ii < N_con) {
        p_con[0][ii] = ii*dcon;
      } else {
        p_con[0][ii] = p_fpt[0][ii-N_con+1] + tnd_a;
      }
    }

    for (ii = 1; ii <= N_cut+N_con-1; ii++) {

      ll = N_cut;
      mm = N_cut;

      for (jj = 1; jj < N_dtau; jj++){

        tau0 = tnd_a + (jj-1.0)*dtau;
        tau1 = tnd_a + jj*dtau;
        tmtau0 = p_con[0][ii] - tau0;
        tmtau1 = p_con[0][ii] - tau1;

        if (tmtau0 > 0.0) {
          for (kk = ll; kk > 0; kk--) {
            if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
              p_tnd = tnd_uni(tau0, tnd_a, tnd_b);
              p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
              p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
              ll = kk;
              break;
            }
          }
        }

        if (tmtau1 > 0.0) {
          for (kk = mm; kk > 0; kk--) {
            if ((tmtau1 >= p_fpt[0][kk-1]) && (tmtau1 < p_fpt[0][kk])) {
              p_tnd = tnd_uni(tau1, tnd_a, tnd_b);
              p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
              p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
              mm = kk;
              break;
            }
          }
        }

      }
    }
  }

  /* truncated normal non-decision time distribution */
  if (tnd_dist == 2) {

    tnd_a = tnd_mu - tnd_range*tnd_sig;
    tnd_a = std::max(tnd_a, 0.0);
    tnd_b = tnd_mu + tnd_range*tnd_sig;
    tnd_b = std::min(tnd_b, p_fpt[0][N_cut]);
    dtau = (tnd_b - tnd_a)/(N_dtau - 1.0);

    if (tnd_a > 0.0) {

      dcon = (tnd_a)/(N_con - 1.0);
      for (ii = 0; ii <= N_cut + N_con - 1; ii++) {
        if (ii < N_con) {
          p_con[0][ii] = ii*dcon;
        } else {
          p_con[0][ii] = p_fpt[0][ii-N_con+1] + tnd_a;
        }
      }

      for (ii = 1; ii <= N_cut+N_con-1; ii++) {

        ll = N_cut;
        mm = N_cut;

        for (jj = 1; jj < N_dtau; jj++){

          tau0 = tnd_a + (jj-1.0)*dtau;
          tau1 = tnd_a + jj*dtau;
          tmtau0 = p_con[0][ii] - tau0;
          tmtau1 = p_con[0][ii] - tau1;

          if (tmtau0 > 0.0) {
            for (kk = ll; kk > 0; kk--) {
              if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
                p_tnd = tnd_tn(phi, tau0);
                p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
                p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
                ll = kk;
                break;
              }
            }
          }

          if (tmtau1 > 0.0) {
            for (kk = mm; kk > 0; kk--) {
              if ((tmtau1 >= p_fpt[0][kk-1]) && (tmtau1 < p_fpt[0][kk])) {
                p_tnd = tnd_tn(phi, tau1);
                p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
                p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
                mm = kk;
                break;
              }
            }
          }

        }
      }
    } else {

      for (ii = 0; ii <= N_cut; ii++) {
        p_con[0][ii] = p_fpt[0][ii];
      }

      for (ii = 1; ii <= N_cut; ii++) {

        ll = N_cut;
        mm = N_cut;

        for (jj = 1; jj < N_dtau; jj++){

          tau0 = (jj-1.0)*dtau;
          tau1 = jj*dtau;
          tmtau0 = p_con[0][ii] - tau0;
          tmtau1 = p_con[0][ii] - tau1;

          if (tmtau0 > 0.0) {
            for (kk = ll; kk > 0; kk--) {
              if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
                p_tnd = tnd_tn(phi, tau0);
                p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
                p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
                ll = kk;
                break;
              }
            }
          }

          if (tmtau1 > 0.0) {
            for (kk = mm; kk > 0; kk--) {
              if ((tmtau1 >= p_fpt[0][kk-1]) && (tmtau1 < p_fpt[0][kk])) {
                p_tnd = tnd_tn(phi, tau1);
                p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
                p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
                mm = kk;
                break;
              }
            }
          }

        }
      }
    }
  }

  /* inverse gaussian non-decision time distribution */
  if (tnd_dist == 3) {

    for (ii = 0; ii <= N_cut; ii++) {
      p_con[0][ii] = p_fpt[0][ii];
    }

    for (ii = 1; ii <= N_cut; ii++) {

      ll = N_cut;
      mm = N_cut;

      for (jj = 1; jj <= N_cut; jj++){

        tau0 = p_con[0][jj-1];
        tau1 = p_con[0][jj];
        dtau = tau1 - tau0;

        tmtau0 = p_con[0][ii] - tau0;
        tmtau1 = p_con[0][ii] - tau1;

        if (tmtau0 > 0.0) {
          for (kk = ll; kk > 0; kk--) {
            if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
              if (tau0 <= 0.0) {
                p_tnd = 0.0;
              } else {
                p_tnd = tnd_ig(phi, tau0);
              }
              p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
              p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
              ll = kk;
              break;
            }
          }
        }

        if (tmtau1 > 0.0) {
          for (kk = mm; kk > 0; kk--) {
            if ((tmtau1 >= p_fpt[0][kk-1]) && (tmtau1 < p_fpt[0][kk])) {
              if (tau1 <= 0.0) {
                p_tnd = 0.0;
              } else {
                p_tnd = tnd_ig(phi, tau1);
              }
              p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
              p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
              mm = kk;
              break;
            }
          }
        }

      }
    }
  }

  /* ---------------------------------- */
  /* calculate cumulative probabilities */
  /* ---------------------------------- */

  /* include "dummy" time steps into loop count */
  N_print = N_cut;
  if ((tnd_dist == 0) || (tnd_dist == 1)) {
    N_print = N_cut + N_con - 1;
  }
  if ((tnd_dist == 2) && (tnd_a > 0.0)) {
    N_print = N_cut + N_con - 1;
  }

  for (ii = 0; ii <= N_print; ii++) {

    /* include contamination distribution */
    g_prob = 0.5*gg*contamination_probability(phi, p_con[0][ii]);
    p_con[1][ii] = (1.0 - gg)*p_con[1][ii] + g_prob;
    p_con[2][ii] = (1.0 - gg)*p_con[2][ii] + g_prob;

    /* set small values to p_fpt_min */
    if (p_con[1][ii] < p_fpt_min) {
      p_con[1][ii] = p_fpt_min;
    }
    if (p_con[2][ii] < p_fpt_min) {
      p_con[2][ii] = p_fpt_min;
    }

    if (ii == 0) {
      cdf[0][ii] = p_con[1][ii];
      cdf[1][ii] = p_con[2][ii];
    } else {
      cdf[0][ii] = cdf[0][ii-1] + 0.5*(p_con[0][ii] - p_con[0][ii-1])*(p_con[1][ii] + p_con[1][ii-1]);
      cdf[1][ii] = cdf[1][ii-1] + 0.5*(p_con[0][ii] - p_con[0][ii-1])*(p_con[2][ii] + p_con[2][ii-1]);
    }

  }

  /* --------------------------- */
  /* calculate logcdf of rt data */
  /* --------------------------- */

  /* determine size of rt vectors */
  int N_rtu = rtu.size();
  int N_rtl = rtl.size();

  /* initialize loglikelihood variable */
  double logCDF = 0.0;

  /* interpolate to find loglikelihood of data */
  kk = 0;
  for (ii = 0; ii < N_rtu; ii++){
    if ( (rtu[ii] > p_con[0][N_print]) || (rtu[ii] < p_con[0][0]) ) {
      RCDFupp[ii] = cdf[0][N_cut];
      RlogCDFupp[ii] = log(RCDFupp[ii]);
      logCDF += RlogCDFupp[ii];
    } else if ( (rtu[ii] < p_con[0][0]) ) {
      RCDFupp[ii] = p_fpt_min;
      RlogCDFupp[ii] = log(RCDFupp[ii]);
      logCDF += RlogCDFupp[ii];
    } else {
      for (jj = kk; jj < N_print; jj++){
        if ( (rtu[ii] >= p_con[0][jj]) && (rtu[ii] < p_con[0][jj+1]) ) {
          RCDFupp[ii] = cdf[0][jj] + (rtu[ii] - p_con[0][jj]) * (cdf[0][jj+1] - cdf[0][jj]) / (p_con[0][jj+1] - p_con[0][jj]);
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
    if ( (rtl[ii] > p_con[0][N_print]) ) {
      RCDFlow[ii] = cdf[1][N_cut];
      RlogCDFlow[ii] = log(RCDFlow[ii]);
      logCDF += RlogCDFlow[ii];
    } else if ( (rtl[ii] < p_con[0][0]) ) {
      RCDFlow[ii] = p_fpt_min;
      RlogCDFlow[ii] = log(RCDFlow[ii]);
      logCDF += RlogCDFlow[ii];
    } else {
      for (jj = kk; jj < N_print; jj++){
        if ( (rtl[ii] >= p_con[0][jj]) && (rtl[ii] < p_con[0][jj+1]) ) {
          RCDFlow[ii] = cdf[1][jj] + (rtl[ii] - p_con[0][jj]) * (cdf[1][jj+1] - cdf[1][jj]) / (p_con[0][jj+1] - p_con[0][jj]);
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
int Model_7P::rand(double *Rrt, double *phi) const {

  /* initialize loop indicies */
  int ii = 0;

  /* initialize remaining variables */
  double tnd_mean = non_decision_average(phi); /* non-decision mean */
  double tnd_std = non_decision_width(phi); /* non-decision width */
  double tnd = 0.0; /* non-decision time */
  double w_mean = relative_start_average(phi); /* relative start point mean */
  double w_std = relative_start_width(phi); /* relative start point width */
  double ww = 0.0; /* relative start point */
  double v_mean = drift_average(phi); /* mean drift rate */
  double v_std = drift_width(phi); /* width of drift rate distribution */
  double vv = 0.0; /* drift rate */
  double DD = diffusion(phi, 0.0, 0.0); /* diffusion rate */
  double bu = upper_threshold(phi, 0.0); /* upper decision threshold */
  double bl = lower_threshold(phi, 0.0); /* lower decision threshold */
  double sqrtdt = sqrt(dt_); /* calculate square root of time step */
  double tt = 0.0; /* time */
  double xx = 0.0; /* accumulated evidence */
  double uu = 0.0; /* random number value */

  /* seed rng */
  GetRNGstate();

  /* simulate model */
  for (ii = 0; ii < N; ii++) {

    /* reset time */
    tt = 0.0;

    /* generate non-decision time */
    if (tnd_dist == 0) {
      tnd = tnd_mean;
    } else if (tnd_dist == 1) {
      tnd = Rf_runif(tnd_mean - tnd_std, tnd_mean + tnd_std);
    } else if (tnd_dist == 2) {
      do {
        tnd = Rf_rnorm(tnd_mean, tnd_std);
      } while (tnd < 0.0); // Ensure non-negative values
    } else if (tnd_dist == 3) {
      uu = norm_rand();
      uu = uu * uu; // Square of the standard normal random variable
      tnd = tnd_mean + (tnd_mean * tnd_mean * uu) / (2.0 * tnd_std)
        - tnd_mean / (2.0 * tnd_std) * sqrt(4.0 * tnd_mean * tnd_std * uu + tnd_mean * tnd_mean * uu * uu);
      uu = unif_rand(); // Uniform random variable for inversion logic
      if (uu > tnd_mean / (tnd_mean + tnd)) {
        tnd = tnd_mean * tnd_mean / tnd; // Apply inverse Gaussian transformation
      }
    }

    /* Generate start point */
    if (w_dist == 0) {
      ww = w_mean;
    } else if (w_dist == 1) {
      ww = Rf_runif(w_min, w_max); // Uniform distribution within [w_min, w_max]
    } else if (w_dist == 2) {
      do {
        ww = Rf_rnorm(w_mean, (w_max - w_min) / 2.0); // Normal distribution
      } while ((ww < w_min) || (ww > w_max)); // Ensure ww is within bounds
    }
    xx = bl + ww * (bu - bl); // Scale start point to the desired range

    /* Generate drift rate */
    if (v_dist == 0) {
      vv = Rf_runif(v_mean - v_std, v_mean + v_std); // Uniform distribution
    } else if (v_dist == 1) {
      vv = Rf_rnorm(v_mean, v_std); // Normal distribution
    }

    while (tt <= tsim_max) {

      /* update time step */
      tt += dt_;

      /* draw random number for simulating diffusion */
      uu = norm_rand();

      /* calculate accumulated evidence */
      xx += dt_*vv + sqrtdt*DD*uu;

      /* check if accumulated evidence has crossed a decision threshold */
      if (xx >= bu) {
        // rt = t_old + (bu - x_old)/(xx - x_old)*(tt - t_old);
        Rrt[ii] = tnd + tt;
        break;
      } else if (xx <= bl) {
        // rt = t_old + (bl - x_old)/(xx - x_old)*(tt - t_old);
        Rrt[ii] = -tnd - tt;
        break;
      }

    }
  }
  PutRNGstate();

  return 0;

}



/* method used to construct grid for PDF */
int Model_7P::grid_pdf(double *Rrt, double *Rpdf_u, double *Rpdf_l, double *phi) const {

  /* initialize loop indicies */
  int ii = 0, jj = 0, kk = 0, ll = 0, mm = 0, nn = 0;

  /* initialize arrays */
  double p_fpt[3][N_dt]; /* first passage time probability */
  double p_atv[2][N_dt]; /* first passage time probability for pdf with across trial variability in drift rate */
  double p_con[3][N_dt]; /* first passage time probability for pdf convolved with non-decision time distribution */
  double p_ini[N_deps_max]; /* initial probability distribution */
  double p_n[N_deps_max]; /* probability of evidence x at time step n */
  double p_ncn[N_deps_max]; /* Crank-Nicolson modified probability of evidence x at time step n */
  double p_np1[N_deps_max]; /* Crank-Nicolson modified probability of evidence x at time step n */
  double AA[N_deps_max], BB[N_deps_max], CC[N_deps_max], DD[N_deps_max], EE[N_deps_max], FF[N_deps_max]; /* Crank-Nicolson matrix arrays */
  double ws[2]; /* start point array */
  int wi[2]; /* index of start point arrays */

  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < N_dt; jj++) {
      p_fpt[ii][jj] = 0.0;
      p_con[ii][jj] = 0.0;
    }
  }

  for (ii = 0; ii < 2; ii++) {
    for (jj = 0; jj < N_dt; jj++) {
      p_atv[ii][jj] = 0.0;
    }
  }

  for (ii = 0; ii < N_deps_max; ii++) {
    p_ini[ii] = 0.0;
    p_n[ii] = 0.0;
    p_ncn[ii] = 0.0;
    p_np1[ii] = 0.0;
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
  double v_avg = drift_average(phi);
  double v_sd = drift_width(phi);
  double vv = 0.0;
  double dv = 0.0;
  double sigma = diffusion(phi, 0.0, 0.0);
  double sigma2 = sigma*sigma;
  double deps = 1.0/(N_deps - 1.0);
  double eps = 0.0;
  double p_ini_int = 0.0;
  double dt_base = 0.0, dt_ini = 0.0, dt = 0.0;
  double t_in = 0.0;
  double bu = upper_threshold(phi, 0.0), bl = lower_threshold(phi, 0.0);
  double ww = 0.0;
  double w_avg = relative_start_average(phi);
  double w_sig = relative_start_width(phi);
  double w_a = 0.0;
  double w_b = 0.0;
  double tt = 0.0;
  double alpha = 0.0;
  double gamma = 0.0;
  double WW = 0.0;
  double s_ini = bu - bl;
  double int_prob = 0.0;
  double gg = contamination_strength(phi);
  double g_prob = 0.0;
  double tnd_mu = non_decision_average(phi);
  double tnd_sig = non_decision_width(phi);
  double tnd_a = 0.0;
  double tnd_b = 0.0;
  double dcon = 0.0;
  double tau0 = 0.0;
  double tau1 = 0.0;
  double dtau = 0.0;
  double tmtau0 = 0.0;
  double tmtau1 = 0.0;
  double p_tnd = 0.0;
  int N_cut = 0;
  int N_print = 0;

  /* calculate drift rate step */
  if (v_dist == 0) {
    dv = ( (v_avg + v_sd) - (v_avg - v_sd) ) / (N_dv - 1.0);
  } else if (v_dist == 1) {
    dv = ( (v_avg + v_range*v_sd) - (v_avg - v_range*v_sd) ) / (N_dv - 1.0);
  }

  /* ------------------------------------------- */
  /* calculate relative start point distribution */
  /* ------------------------------------------- */

  /* delta relative start distribution */
  if (w_dist == 0) {

    ww = w_avg;

    /* calculate probability around start point, then determine nearest array elements to delta distribution */
    wi[0] = int( ceil(ww/deps) );
    wi[1] = int( floor(ww/deps) );

    /* determine distance between delta distribution and nearest array elements */
    ws[0] = fabs( ww - deps*wi[0] );
    ws[1] = fabs( ww - deps*wi[1] );

    /* numerically approximate delta function, distribute between nearest array elements */
    if (wi[0] == wi[1]) {
      p_ini[ wi[0] ] = 1.0/deps;
    } else {
      p_ini[ wi[0] ] = (1.0 - ws[0]/deps)/deps;
      p_ini[ wi[1] ] = (1.0 - ws[1]/deps)/deps;
    }

    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      p_ini_int += 0.5*deps*(p_ini[ii-1] + p_ini[ii]);
    }

    for (ii = 0; ii < N_deps; ii++) {
      p_ini[ii] = (1.0/p_ini_int)*p_ini[ii];
    }

  }

  /* uniform relative start distribution */
  if (w_dist == 1) {

    w_a = std::max(w_avg - w_sig, w_min);
    w_b = std::min(w_avg + w_sig, w_max);

    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      if ((w_a >= eps - deps) && (w_a < eps)) {
        if (fabs(w_a - (eps - deps)) <= fabs(eps - w_a)) {
          w_a = eps - deps;
        } else {
          w_a = eps;
        }
        break;
      }
    }

    for (ii = N_deps; ii > 0; ii--) {
      eps = ii*deps;
      if ((w_b >= eps - deps) && (w_b < eps)) {
        if (fabs(w_b - (eps - deps)) <= fabs(eps - w_b)) {
          w_b = eps - deps;
        } else {
          w_b = eps;
        }
        break;
      }
    }

    p_ini[0] = w_uni(eps, w_a, w_b);
    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      p_ini[ii] = w_uni(eps, w_a, w_b);
      p_ini_int += 0.5*deps*(p_ini[ii-1] + p_ini[ii]);
    }

    for (ii = 0; ii < N_deps; ii++) {
      p_ini[ii] = 1.0/p_ini_int*p_ini[ii];
    }

  }

  /* truncated normal relative start distribution */
  if (w_dist == 2) {

    w_a = w_min;
    w_b = w_max;

    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      if ((w_a >= eps - deps) && (w_a < eps)) {
        if (fabs(w_a - (eps - deps)) <= fabs(eps - w_a)) {
          w_a = eps - deps;
        } else {
          w_a = eps;
        }
        break;
      }
    }

    for (ii = N_deps; ii > 0; ii--) {
      eps = ii*deps;
      if ((w_b >= eps - deps) && (w_b < eps)) {
        if (fabs(w_b - (eps - deps)) <= fabs(eps - w_b)) {
          w_b = eps - deps;
        } else {
          w_b = eps;
        }
        break;
      }
    }

    for (ii = 1; ii < N_deps; ii++) {
      eps = ii*deps;
      if ((eps >= w_a) && (eps <= w_b)) {
        p_ini[ii] = w_tn(phi, eps, w_min, w_max);
      }
      p_ini_int += 0.5*deps*(p_ini[ii-1] + p_ini[ii]);
    }

    for (ii = 0; ii < N_deps; ii++) {
      p_ini[ii] = (1.0/p_ini_int)*p_ini[ii];
    }

  }

  /* -------------------------------------------------------------------------------------------- */
  /* calculate first passage time distribution N_dv times for drift rate across trial variability */
  /* -------------------------------------------------------------------------------------------- */

  /* calculate time step sizes */
  dt_base = approx_dt(phi, dt_scale);
  if (dt_base > dt_max) {
    dt_base = dt_max;
  }
  dt_ini = dt_ini_scale*dt_base;
  dt = dt_ini;

  /* set non-decision width to dt_base if smaller */
  tnd_sig = std::max(tnd_sig, dt_base);
  tnd_a = tnd_mu - tnd_sig;
  tnd_b = tnd_mu + tnd_sig;

  for (nn = 0; nn < N_dv; nn++) {

    /* reset arrays */
    for (ii = 1; ii < 3; ii++) {
      for (jj = 0; jj < N_dt; jj++) {
        p_fpt[ii][jj] = 0.0;
      }
    }

    for (ii = 0; ii < N_deps_max; ii++) {
      p_n[ii] = p_ini[ii];
      p_ncn[ii] = 0.0;
      p_np1[ii] = 0.0;
      AA[ii] = 0.0;
      BB[ii] = 0.0;
      CC[ii] = 0.0;
      DD[ii] = 0.0;
      EE[ii] = 0.0;
      FF[ii] = 0.0;
    }

    /* reset constants */
    if (v_dist == 0){
      vv = v_avg - v_sd + nn*dv;
    } else if (v_dist == 1) {
      vv = v_avg - v_range*v_sd + nn*dv;
    }

    tt = 0.0;
    alpha = 0.0;
    gamma = 0.0;
    WW = 0.0;
    int_prob = 0.0;
    g_prob = 0.0;

    /* calculate first passage time distribution */
    for (ii = 0; ii < N_dt - 1; ii++) {

      /* set time step */
      if (ii < N_dt_ini) {
        dt = dt_ini;
      } else if ((ii >= N_dt_ini) && (ii < N_dt_ini + N_dt_scaleup)) {
        dt = dt_ini + (ii + 1.0 - N_dt_ini)*(dt_base - dt_ini)/(N_dt_scaleup);
      } else {
        dt = dt_base;
      }
      dt = modify_dt(phi, tt)*dt;
      tt += dt;

      /* invert matrix using tridiagonal matrix algorithm */
      p_np1[0] = 0.0;
      p_np1[1] = 0.0;

      for (jj = 1; jj < N_deps-1; jj++) {

        p_np1[jj+1] = 0.0;

        alpha = dt/(4.0*s_ini*deps);
        gamma = alpha/(s_ini*deps);

        AA[jj] = -alpha*vv - gamma*sigma2;
        BB[jj] = 1.0 + 2.0*gamma*sigma2;
        CC[jj] = alpha*vv - gamma*sigma2;

        DD[jj] = alpha*vv + gamma*sigma2;
        EE[jj] = 1.0 - 2.0*gamma*sigma2;
        FF[jj] = -alpha*vv + gamma*sigma2;

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
      p_fpt[0][ii+1] = tt;
      p_fpt[1][ii+1] = fabs( 0.5 * sigma2 * (4.0*p_np1[N_deps-2] - p_np1[N_deps-3])/(2.0*deps*s_ini*s_ini) );
      p_fpt[2][ii+1] = fabs( 0.5 * sigma2 * (4.0*p_np1[1] - p_np1[2])/(2.0*deps*s_ini*s_ini) );

      /* determine integrated probability */
      int_prob += p_fpt[1][ii+1]*dt + p_fpt[2][ii+1]*dt;

      if (v_dist == 0) {
        if ( (nn == 0) || (nn == N_dv - 1) ) {
          p_atv[0][ii+1] += dv*v_uni(phi, vv)*p_fpt[1][ii+1]/2.0;
          p_atv[1][ii+1] += dv*v_uni(phi, vv)*p_fpt[2][ii+1]/2.0;
        } else {
          p_atv[0][ii+1] += dv*v_uni(phi, vv)*p_fpt[1][ii+1];
          p_atv[1][ii+1] += dv*v_uni(phi, vv)*p_fpt[2][ii+1];
        }
      } else if (v_dist == 1) {
        if ( (nn == 0) || (nn == N_dv - 1) ) {
          p_atv[0][ii+1] += dv*v_norm(phi, vv)*p_fpt[1][ii+1]/2.0;
          p_atv[1][ii+1] += dv*v_norm(phi, vv)*p_fpt[2][ii+1]/2.0;
        } else {
          p_atv[0][ii+1] += dv*v_norm(phi, vv)*p_fpt[1][ii+1];
          p_atv[1][ii+1] += dv*v_norm(phi, vv)*p_fpt[2][ii+1];
        }
      }

      /* determine index of highest time point */
      N_cut = std::max(ii+1, N_cut);

      /* stop solver if time above rt_max or probability below p_fpt_min */
      if (p_fpt[0][ii+1] > rt_max) {
        break;
      } else if ((int_prob > int_prob_min) && (p_fpt[1][ii+1] <= p_fpt_min) && (p_fpt[2][ii+1] <= p_fpt_min)) {
        break;
      }

    }

  }

  /* ------------------------------------------ */
  /* convolve pdf and non-decision distribution */
  /* ------------------------------------------ */

  /* delta non-decision time distribution */
  if (tnd_dist == 0) {

    dcon = tnd_mu/(N_con - 1.0);
    for (ii = 0; ii <= N_cut + N_con - 1; ii++) {
      if (ii < N_con) {
        p_con[0][ii] = ii*dcon;
      } else {
        p_con[0][ii] = p_fpt[0][ii-N_con+1] + tnd_mu;
      }
    }

    ll = 1;

    for (ii = 1; ii <= N_cut+N_con-1; ii++) {

      tmtau0 = p_con[0][ii] - tnd_mu;

      if (tmtau0 > 0.0) {
        for (kk = ll; kk <= N_cut; kk++) {
          if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
            p_con[1][ii] = interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
            p_con[2][ii] = interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
            ll = kk;
            break;
          }
        }
      }

    }
  }

  /* uniform non-decision time distribution */
  if (tnd_dist == 1) {

    tnd_b = std::min(tnd_b, p_fpt[0][N_cut]);
    dtau = (tnd_b - tnd_a)/(N_dtau - 1.0);

    dcon = tnd_a/(N_con - 1.0);
    for (ii = 0; ii <= N_cut + N_con - 1; ii++) {
      if (ii < N_con) {
        p_con[0][ii] = ii*dcon;
      } else {
        p_con[0][ii] = p_fpt[0][ii-N_con+1] + tnd_a;
      }
    }

    for (ii = 1; ii <= N_cut+N_con-1; ii++) {

      ll = N_cut;
      mm = N_cut;

      for (jj = 1; jj < N_dtau; jj++){

        tau0 = tnd_a + (jj-1.0)*dtau;
        tau1 = tnd_a + jj*dtau;
        tmtau0 = p_con[0][ii] - tau0;
        tmtau1 = p_con[0][ii] - tau1;

        if (tmtau0 > 0.0) {
          for (kk = ll; kk > 0; kk--) {
            if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
              p_tnd = tnd_uni(tau0, tnd_a, tnd_b);
              p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
              p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
              ll = kk;
              break;
            }
          }
        }

        if (tmtau1 > 0.0) {
          for (kk = mm; kk > 0; kk--) {
            if ((tmtau1 >= p_fpt[0][kk-1]) && (tmtau1 < p_fpt[0][kk])) {
              p_tnd = tnd_uni(tau1, tnd_a, tnd_b);
              p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
              p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
              mm = kk;
              break;
            }
          }
        }

      }
    }
  }

  /* truncated normal non-decision time distribution */
  if (tnd_dist == 2) {

    tnd_a = tnd_mu - tnd_range*tnd_sig;
    tnd_a = std::max(tnd_a, 0.0);
    tnd_b = tnd_mu + tnd_range*tnd_sig;
    tnd_b = std::min(tnd_b, p_fpt[0][N_cut]);
    dtau = (tnd_b - tnd_a)/(N_dtau - 1.0);

    if (tnd_a > 0.0) {

      dcon = (tnd_a)/(N_con - 1.0);
      for (ii = 0; ii <= N_cut + N_con - 1; ii++) {
        if (ii < N_con) {
          p_con[0][ii] = ii*dcon;
        } else {
          p_con[0][ii] = p_fpt[0][ii-N_con+1] + tnd_a;
        }
      }

      for (ii = 1; ii <= N_cut+N_con-1; ii++) {

        ll = N_cut;
        mm = N_cut;

        for (jj = 1; jj < N_dtau; jj++){

          tau0 = tnd_a + (jj-1.0)*dtau;
          tau1 = tnd_a + jj*dtau;
          tmtau0 = p_con[0][ii] - tau0;
          tmtau1 = p_con[0][ii] - tau1;

          if (tmtau0 > 0.0) {
            for (kk = ll; kk > 0; kk--) {
              if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
                p_tnd = tnd_tn(phi, tau0);
                p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
                p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
                ll = kk;
                break;
              }
            }
          }

          if (tmtau1 > 0.0) {
            for (kk = mm; kk > 0; kk--) {
              if ((tmtau1 >= p_fpt[0][kk-1]) && (tmtau1 < p_fpt[0][kk])) {
                p_tnd = tnd_tn(phi, tau1);
                p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
                p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
                mm = kk;
                break;
              }
            }
          }

        }
      }
    } else {

      for (ii = 0; ii <= N_cut; ii++) {
        p_con[0][ii] = p_fpt[0][ii];
      }

      for (ii = 1; ii <= N_cut; ii++) {

        ll = N_cut;
        mm = N_cut;

        for (jj = 1; jj < N_dtau; jj++){

          tau0 = (jj-1.0)*dtau;
          tau1 = jj*dtau;
          tmtau0 = p_con[0][ii] - tau0;
          tmtau1 = p_con[0][ii] - tau1;

          if (tmtau0 > 0.0) {
            for (kk = ll; kk > 0; kk--) {
              if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
                p_tnd = tnd_tn(phi, tau0);
                p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
                p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
                ll = kk;
                break;
              }
            }
          }

          if (tmtau1 > 0.0) {
            for (kk = mm; kk > 0; kk--) {
              if ((tmtau1 >= p_fpt[0][kk-1]) && (tmtau1 < p_fpt[0][kk])) {
                p_tnd = tnd_tn(phi, tau1);
                p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
                p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
                mm = kk;
                break;
              }
            }
          }

        }
      }
    }
  }

  /* inverse gaussian non-decision time distribution */
  if (tnd_dist == 3) {

    for (ii = 0; ii <= N_cut; ii++) {
      p_con[0][ii] = p_fpt[0][ii];
    }

    for (ii = 1; ii <= N_cut; ii++) {

      ll = N_cut;
      mm = N_cut;

      for (jj = 1; jj <= N_cut; jj++){

        tau0 = p_con[0][jj-1];
        tau1 = p_con[0][jj];
        dtau = tau1 - tau0;

        tmtau0 = p_con[0][ii] - tau0;
        tmtau1 = p_con[0][ii] - tau1;

        if (tmtau0 > 0.0) {
          for (kk = ll; kk > 0; kk--) {
            if ((tmtau0 >= p_fpt[0][kk-1]) && (tmtau0 < p_fpt[0][kk])) {
              if (tau0 <= 0.0) {
                p_tnd = 0.0;
              } else {
                p_tnd = tnd_ig(phi, tau0);
              }
              p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
              p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau0, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
              ll = kk;
              break;
            }
          }
        }

        if (tmtau1 > 0.0) {
          for (kk = mm; kk > 0; kk--) {
            if ((tmtau1 >= p_fpt[0][kk-1]) && (tmtau1 < p_fpt[0][kk])) {
              if (tau1 <= 0.0) {
                p_tnd = 0.0;
              } else {
                p_tnd = tnd_ig(phi, tau1);
              }
              p_con[1][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[0][kk-1], p_atv[0][kk]);
              p_con[2][ii] += 0.5*dtau*p_tnd*interp(tmtau1, p_fpt[0][kk-1], p_fpt[0][kk], p_atv[1][kk-1], p_atv[1][kk]);
              mm = kk;
              break;
            }
          }
        }

      }
    }
  }

  /* ------------------------- */
  /* output probabilities to R */
  /* ------------------------- */

  /* include "dummy" time steps into loop count */
  N_print = N_cut;
  if ((tnd_dist == 0) || (tnd_dist == 1)) {
    N_print = N_cut + N_con - 1;
  }
  if ((tnd_dist == 2) && (tnd_a > 0.0)) {
    N_print = N_cut + N_con - 1;
  }

  for (ii = 0; ii <= N_print; ii++) {

    /* include contamination distribution */
    g_prob = 0.5*gg*contamination_probability(phi, p_con[0][ii]);
    p_con[1][ii] = (1.0 - gg)*p_con[1][ii] + g_prob;
    p_con[2][ii] = (1.0 - gg)*p_con[2][ii] + g_prob;

    /* set small values to p_fpt_min */
    if (p_con[1][ii] < p_fpt_min) {
      p_con[1][ii] = p_fpt_min;
    }
    if (p_con[2][ii] < p_fpt_min) {
      p_con[2][ii] = p_fpt_min;
    }

    Rrt[ii+1] = p_con[0][ii];
    Rpdf_u[ii+1] = p_con[1][ii];
    Rpdf_l[ii+1] = p_con[2][ii];

  }

  return 0;

}
