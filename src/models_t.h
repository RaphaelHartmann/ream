/* file: model_functions.cpp
 Functions for defining model.
 Author: Mathew Murrow and Raphael Hartmann
 Date: Sep 02, 2024 */

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

#ifndef MODELS_T_H
#define MODELS_T_H

#include "Model_T.h"



class DMC : public Model_T {
protected:

  /* function for the non-decision time */
  double non_decision(const double phi[12]) const override {
    return phi[0];
  }

  /* function for the start point */
  double relative_start(const double phi[12]) const override {
    return phi[1];
  }

  /* function for the drift rate */
  double drift(const double phi[12], double t) const override {
    double e = exp(1.0); /* natural exponent */
    double t_floor = 1.0e-9; /* min t to prevent singularity at zero */
    double con = phi[2]; /* congruency (1.0 is congruent, -1.0 is incongruent) */
    double A = phi[3]; /* peak amplitude of automatic process drift rate */
    double tau = phi[4]; /* characteristic time of automatic process drift rate */
    double alpha = phi[5]; /* shape of automatic process drift rate */
    double mu_c = phi[6]; /* drift rate for controlled process */
    double v = con*A*exp(-t/tau) * pow( (t*e)/( (alpha-1.0)*tau )  , alpha-1.0 ) * ( (alpha-1.0)/(t + t_floor) - 1.0/tau ) + mu_c; /* total DMC drift rate */
    return v;
  }

  /* function for the diffusion rate */
  double diffusion(const double phi[12], double x, double t) const override {
    return phi[7];
  }

  /* function for the upper threshold */
  double upper_threshold(const double phi[12], double t) const override {
    return phi[8];
  }

  /* function for the lower threshold */
  double lower_threshold(const double phi[12], double t) const override {
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

};



class PAM : public Model_T {
protected:

  /* method for the non-decision time of process 1 */
  double non_decision(const double phi[11]) const override {
    return phi[0];
  }

  /* method for the start point of process 1 */
  double relative_start(const double phi[11]) const override {
    return phi[1];
  }

  /* method for the drift rate of process 1 */
  double drift(const double phi[11], double t) const override {
    double a_inner = 1.0;
    double a_outer = 1.0;
    double a_target = 1.0;
    double p_outer = phi[2];
    double p_inner = phi[3];
    double p_target = phi[4];
    double t_s = phi[5];
    double v = 0.0;

    if (t >= t_s) {
      a_inner = 0.0;
      a_outer = 0.0;
    }

    v = 2.0*a_outer*p_outer + 2.0*a_inner*p_inner + a_target*p_target;

    return v;
  }

  /* method for the diffusion rate of process 1 */
  double diffusion(const double phi[11], double x, double t) const override {
    return phi[6];
  }

  /* method for the upper threshold of process 1 */
  double upper_threshold(const double phi[11], double t) const override {
    return phi[7];
  }

  /* method for the lower threshold of process 1 */
  double lower_threshold(const double phi[11], double t) const override {
    return -phi[7];
  }

  /* method for the contamination strength of process 1 */
  double contamination_strength(const double phi[11]) const override {
    return phi[8];
  }

  /* method for the contamination probability distribution of process 1 */
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
    double t_s = phi[5];
    double range = 0.1;
    double dt_mod = 1.0;

    if ((t >= t_s - range) && (t <= t_s + range)) {
      dt_mod = 0.1;
    }

    return dt_mod;
  }

};



// class DSTP : public Model_T {
// protected:
//
//   /* method for the non-decision time of process 1 */
//   double non_decision(const double phi[15]) const override {
//     return phi[0];
//   }
//
//   /* method for the start point of process 1 */
//   double relative_start(const double phi[15]) const override {
//     return phi[1];
//   }
//
//   /* method for the drift rate of process 1 */
//   double drift(const double phi[15], double t) const override {
//   double w2 = phi[2]; /* relative start point for process 2 */
//   double v2 = phi[3]; /* drift rate for process 2 */
//   double sigma2 = phi[4]; /* diffusion rate for process 2 */
//   double a2 = phi[5]; /* threshold separation for process 2 */
//   double z2 = w2*a2; /* start point for process 2 */
//   double mu_tar = phi[6]; /* drift rate for target in process 1 */
//   double mu_fl = phi[7]; /* drift rate for flanker in process 1 */
//   double mu_rs2 = phi[8]; /* stage 2 drift rate for process 1 */
//   double v = 0.0; /* total summed drift rate for process 1 */
//   int kk = 0, N_k = phi[9];
//
//   /* calculate probability p of process 2 crossing upper and lower threhsolds */
//   double p_upper = ( exp(-2.0*(-v2)*a2/(sigma2*sigma2)) - exp(-2.0*(-v2)*(a2-z2)/(sigma2*sigma2)) ) / ( exp(-2.0*(-v2)*a2/(sigma2*sigma2)) - 1.0 );
//   double p_lower = ( exp(-2.0*v2*a2/(sigma2*sigma2)) - exp(-2.0*v2*z2/(sigma2*sigma2)) ) / ( exp(-2.0*v2*a2/(sigma2*sigma2)) - 1.0 );
//
//   /* calculate cumulative probability g for upper and lower threhsolds */
//   double g_upper = 0.0;
//   double g_lower = 0.0;
//   double weight = 0.0;
//
//   for (kk = 1; kk < N_k + 1; kk++){
//     g_upper += 2.0*kk*sin(kk*M_PI*(a2-z2)/a2)*exp(-0.5*t*((-v2)*(-v2)/(sigma2*sigma2) + (M_PI*M_PI)*(kk*kk)*(sigma2*sigma2)/(a2*a2))) / ((-v2)*(-v2)/(sigma2*sigma2) + (M_PI*M_PI)*(kk*kk)*(sigma2*sigma2)/(a2*a2));
//     g_lower += 2.0*kk*sin(kk*M_PI*z2/a2)*exp(-0.5*t*((v2*v2)/(sigma2*sigma2) + (M_PI*M_PI)*(kk*kk)*(sigma2*sigma2)/(a2*a2))) / ((v2*v2)/(sigma2*sigma2) + (M_PI*M_PI)*(kk*kk)*(sigma2*sigma2)/(a2*a2));
//   }
//   g_upper = p_upper - M_PI*(sigma2*sigma2)/(a2*a2)*exp(v2*(a2-z2)/(sigma2*sigma2))*g_upper;
//   g_lower = p_lower - M_PI*(sigma2*sigma2)/(a2*a2)*exp(-v2*z2/(sigma2*sigma2))*g_lower;
//
//   /* calculate weighted drift rate*/
//   weight = g_upper + g_lower;
//   v = (1.0 - weight)*(mu_tar + mu_fl) + weight*mu_rs2;
//
//   return v;
//   }
//
//   /* method for the diffusion rate of process 1 */
//   double diffusion(const double phi[15], double x, double t) const override {
//     return phi[10];
//   }
//
//   /* method for the upper threshold of process 1 */
//   double upper_threshold(const double phi[15], double t) const override {
//     return phi[11];
//   }
//
//   /* method for the lower threshold of process 1 */
//   double lower_threshold(const double phi[15], double t) const override {
//     return -phi[11];
//   }
//
//   /* method for the contamination strength of process 1 */
//   double contamination_strength(const double phi[15]) const override {
//     return phi[12];
//   }
//
//   /* method for the contamination probability distribution of process 1 */
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
//   /* method for locally modifying the time step size */
//   double modify_dt(const double phi[15], double t) const override {
//     return 1.0;
//   }
//
// };



class ETM : public Model_T {
protected:

  /* method for the non-decision time of process 1 */
  double non_decision(const double phi[9]) const override {
    return phi[0];
  }

  /* method for the start point of process 1 */
  double relative_start(const double phi[9]) const override {
    return phi[1];
  }

  /* method for the drift rate of process 1 */
  double drift(const double phi[9], double t) const override {
    return phi[2];
  }

  /* method for the diffusion rate of process 1 */
  double diffusion(const double phi[9], double x, double t) const override {
    return phi[3];
  }

  /* method for the upper threshold of process 1 */
  double upper_threshold(const double phi[9], double t) const override {
    double b = phi[4];
    double tau = pow(10.0, phi[5]);
    double thres = b*exp(-t/tau);
    return thres;
  }

  /* method for the lower threshold of process 1 */
  double lower_threshold(const double phi[9], double t) const override {
    return -upper_threshold(phi, t);
  }

  /* method for the contamination strength of process 1 */
  double contamination_strength(const double phi[9]) const override {
    return phi[6];
  }

  /* method for the contamination probability distribution of process 1 */
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



class LTM : public Model_T {
protected:

  /* method for the non-decision time of process 1 */
  double non_decision(const double phi[9]) const override {
    return phi[0];
  }

  /* method for the start point of process 1 */
  double relative_start(const double phi[9]) const override {
    return phi[1];
  }

  /* method for the drift rate of process 1 */
  double drift(const double phi[9], double t) const override {
    return phi[2];
  }

  /* method for the diffusion rate of process 1 */
  double diffusion(const double phi[9], double x, double t) const override {
    return phi[3];
  }

  /* method for the upper threshold of process 1 */
  double upper_threshold(const double phi[9], double t) const override {
    double b = phi[4];
    double m = phi[5];
    double thres = b - m*t;
    return thres;
  }

  /* method for the lower threshold of process 1 */
  double lower_threshold(const double phi[9], double t) const override {
    return -upper_threshold(phi, t);
  }

  /* method for the contamination strength of process 1 */
  double contamination_strength(const double phi[9]) const override {
    return phi[6];
  }

  /* method for the contamination probability distribution of process 1 */
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



class RDMC : public Model_T {
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
  double drift(const double phi[11], double t) const override {
    double A_0 = phi[2];
    double k = phi[3];
    double d_a = phi[4];
    double d_c = phi[5];
    double w_a = A_0*exp(-k*t);
    double w_c = 1.0 - w_a;
    double v = w_a*d_a + w_c*d_c;
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
    return 1.0;
  }

};



class RTM : public Model_T {
protected:

  /* method for the non-decision time of process 1 */
  double non_decision(const double phi[10]) const override {
    return phi[0];
  }

  /* method for the start point of process 1 */
  double relative_start(const double phi[10]) const override {
    return phi[1];
  }

  /* method for the drift rate of process 1 */
  double drift(const double phi[10], double t) const override {
    return phi[2];
  }

  /* method for the diffusion rate of process 1 */
  double diffusion(const double phi[10], double x, double t) const override {
    return phi[3];
  }

  /* method for the upper threshold of process 1 */
  double upper_threshold(const double phi[10], double t) const override {
    double b0 = phi[4];
    double kappa = phi[5];
    double t0 = phi[6];
    double thres = 0.5*b0*(1 - kappa*t/(t+t0));
    return thres;
  }

  /* method for the lower threshold of process 1 */
  double lower_threshold(const double phi[10], double t) const override {
    return -upper_threshold(phi, t);
  }

  /* method for the contamination strength of process 1 */
  double contamination_strength(const double phi[10]) const override {
    return phi[7];
  }

  /* method for the contamination probability distribution of process 1 */
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



class SDDM : public Model_T {
protected:

  /* method for the non-decision time of process 1 */
  double non_decision(const double phi[8]) const override {
    return phi[0];
  }

  /* method for the start point of process 1 */
  double relative_start(const double phi[8]) const override {
    return phi[1];
  }

  /* method for the drift rate of process 1 */
  double drift(const double phi[8], double t) const override {
    return phi[2];
  }

  /* method for the diffusion rate of process 1 */
  double diffusion(const double phi[8], double x, double t) const override {
    return phi[3];
  }

  /* method for the upper threshold of process 1 */
  double upper_threshold(const double phi[8], double t) const override {
    return phi[4];
  }

  /* method for the lower threshold of process 1 */
  double lower_threshold(const double phi[8], double t) const override {
    return -phi[4];
  }

  /* method for the contamination strength of process 1 */
  double contamination_strength(const double phi[8]) const override {
    return phi[5];
  }

  /* method for the contamination probability distribution of process 1 */
  double contamination_probability(const double phi[8], double t) const override {
    double gl = phi[6];
    double gu = phi[7];
    double pg = 0.0;
    if ((t >= gl) && (t <= gu)) {
      pg = 1.0/(gu - gl);
    }
    return pg;
  }

  /* method for locally modifying the time step size */
  double modify_dt(const double phi[8], double t) const override {
    return 1.0;
  }

};



class SSP : public Model_T {
protected:

  /* method for the non-decision time of process 1 */
  double non_decision(const double phi[16]) const override {
    return phi[0];
  }

  /* method for the start point of process 1 */
  double relative_start(const double phi[16]) const override {
    return phi[1];
  }

  /* method for the drift rate of process 1 */
  double drift(const double phi[16], double t) const override {
    double sd_a0 = phi[2];
    double r_d = phi[3];
    double sd_a = sd_a0 - r_d * t;
    double c = phi[4];
    if (sd_a < 0.001) {
      sd_a = 0.001;
    }
    double lb_target = phi[5];
    double ub_target = phi[6];
    double lb_inner = ub_target;
    double ub_inner = phi[7];
    double lb_outer = ub_inner;
    double p_target = phi[8];
    double p_inner = c*phi[9];
    double p_outer = c*phi[10];
    double v_target = ncdf(ub_target/sd_a) - ncdf(lb_target/sd_a);
    double v_inner = 2.0*(ncdf(ub_inner/sd_a) - ncdf(lb_inner/sd_a));
    double v_outer = 2.0*(1.0 - ncdf(lb_outer/sd_a));
    double v = p_target*v_target + p_inner*v_inner + p_outer*v_outer;
    return v;
  }

  /* method for the diffusion rate of process 1 */
  double diffusion(const double phi[16], double x, double t) const override {
    return phi[11];
  }

  /* method for the upper threshold of process 1 */
  double upper_threshold(const double phi[16], double t) const override {
    return phi[12];
  }

  /* method for the lower threshold of process 1 */
  double lower_threshold(const double phi[16], double t) const override {
    return -phi[12];
  }

  /* method for the contamination strength of process 1 */
  double contamination_strength(const double phi[16]) const override {
    return phi[13];
  }

  /* method for the contamination probability distribution of process 1 */
  double contamination_probability(const double phi[16], double t) const override {
    double gl = phi[14];
    double gu = phi[15];
    double pg = 0.0;
    if ((t >= gl) && (t <= gu)) {
      pg = 1.0/(gu - gl);
    }
    return pg;
  }

  /* method for locally modifying the time step size */
  double modify_dt(const double phi[16], double t) const override {
    return 1.0;
  }

};



class WTM : public Model_T {
protected:

  /* method for the non-decision time of process 1 */
  double non_decision(const double phi[11]) const override {
    return phi[0];
  }

  /* method for the start point of process 1 */
  double relative_start(const double phi[11]) const override {
    return phi[1];
  }

  /* method for the drift rate of process 1 */
  double drift(const double phi[11], double t) const override {
    return phi[2];
  }

  /* method for the diffusion rate of process 1 */
  double diffusion(const double phi[11], double x, double t) const override {
    return phi[3];
  }

  /* method for the upper threshold of process 1 */
  double upper_threshold(const double phi[11], double t) const override {
    double b = phi[4];
    double lamb = pow(10.0, phi[5]);
    double kappa = pow(10.0, phi[6]);
    double c = phi[7];
    double thres = b - 0.5*b*(1.0 - c)*( 1.0 - exp( -pow( t/lamb, kappa) ) );
    return thres;
  }

  /* method for the lower threshold of process 1 */
  double lower_threshold(const double phi[11], double t) const override {
    return -upper_threshold(phi, t);
  }

  /* method for the contamination strength of process 1 */
  double contamination_strength(const double phi[11]) const override {
    return phi[8];
  }

  /* method for the contamination probability distribution of process 1 */
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
    return 1.0;
  }

};



class CSTM_T : public Model_T {
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
  double drift(const double phi[100], double t) const override {
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
