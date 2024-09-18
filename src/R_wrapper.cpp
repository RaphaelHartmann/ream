/* file: model_functions.cpp
 Functions for defining model.
 Author: Raphael Hartmann and Mathew Murrow
 Date: Sep 02, 2024 */

/* -------------------------------------------------- */
/* -------------------------------------------------- */
/* -------------------------------------------------- */

#define R_NO_REMAP
// #include <chrono>
// #include <thread>
#include "Model.h"
#include "models_t.h"
#include "models_tx.h"
#include "models_tw.h"
#include "tools.h"
#include <R.h>
#include <Rinternals.h>


/* global variables */
const char *PHI;
const char *ModelName;


  // for RAND
int N;
double dt_;

  // for PDF and CDF
const char *RTL;
const char *RTU;
int N_deps;
double dt_scale;
double rt_max;
int N_rtl;
int N_rtu;
int N_phi;

  // for likelihood
// const char *OUTPUT2;
// const char *OUTPUT3;




/* Definition of the createModel function */
std::unique_ptr<Model> createModel(const char* modelName) {
  std::string modelNameStr(modelName);

  if (modelNameStr == "DMC") {
    return std::make_unique<DMC>();
  } else if (modelNameStr == "CDSTP") {
    return std::make_unique<CDSTP>();
  } else if (modelNameStr == "ETM") {
    return std::make_unique<ETM>();
  } else if (modelNameStr == "LTM") {
    return std::make_unique<LTM>();
  } else if (modelNameStr == "PAM") {
    return std::make_unique<PAM>();
  } else if (modelNameStr == "RDMC") {
    return std::make_unique<RDMC>();
  } else if (modelNameStr == "RTM") {
    return std::make_unique<RTM>();
  } else if (modelNameStr == "SDDM") {
    return std::make_unique<SDDM>();
  } else if (modelNameStr == "SDPM") {
    return std::make_unique<SDPM>();
  } else if (modelNameStr == "SSP") {
    return std::make_unique<SSP>();
  } else if (modelNameStr == "UGM") {
    return std::make_unique<UGM>();
  } else if (modelNameStr == "WTM") {
    return std::make_unique<WTM>();
  } else if (modelNameStr == "LIMF") {
    return std::make_unique<LIMF>();
  } else if (modelNameStr == "LIM") {
    return std::make_unique<LIM>();
  } else if (modelNameStr == "UGMF") {
    return std::make_unique<UGMF>();
  } else if (modelNameStr == "WDSTP") {
    return std::make_unique<WDSTP>();
  } else if (modelNameStr == "CSTM_T") {
    return std::make_unique<CSTM_T>();
  } else if (modelNameStr == "CSTM_TX") {
    return std::make_unique<CSTM_TX>();
  } else if (modelNameStr == "CSTM_TW") {
    return std::make_unique<CSTM_TW>();
  } else {
    Rprintf("unknown model name");
    return nullptr;
  }

}




extern "C" {

	SEXP PDF(SEXP re, SEXP in, SEXP re_l, SEXP re_u, SEXP ch) {

	  /* define input variables */

	  ModelName = R_CHAR(STRING_ELT(ch, 0));

	  N_deps = INTEGER(in)[0];
	  N_rtl = INTEGER(in)[1];
	  N_rtu = INTEGER(in)[2];
	  dt_scale = REAL(re)[0];
	  rt_max = REAL(re)[1];

	  std::vector<double> rtl(N_rtl);
	  std::vector<double> rtu(N_rtu);
	  for (int i = 0; i < N_rtl; i++) {
	    rtl[i] = REAL(re_l)[i];
	  }
	  for (int i = 0; i < N_rtu; i++) {
	    rtu[i] = REAL(re_u)[i];
	  }
	  N_phi = INTEGER(in)[3];
	  double *phi = (double*)R_Calloc(N_phi, double);
	  for(int i=0; i<N_phi; i++) {
	    phi[i] = REAL(re)[i+2];
	  }


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
	  SEXP likl = PROTECT(Rf_allocVector(REALSXP, N_rtl));
	  outCnt++;
	  SEXP liku = PROTECT(Rf_allocVector(REALSXP, N_rtu));
	  outCnt++;
	  SEXP llikl = PROTECT(Rf_allocVector(REALSXP, N_rtl));
	  outCnt++;
	  SEXP lliku = PROTECT(Rf_allocVector(REALSXP, N_rtu));
	  outCnt++;
	  SEXP llsum = PROTECT(Rf_allocVector(REALSXP, 1));
	  outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rlikl = REAL(likl);
		double *Rliku = REAL(liku);
		double *Rllikl = REAL(llikl);
		double *Rlliku = REAL(lliku);
		double *Rllsum = REAL(llsum);


		/* model creation */
		auto model = createModel(ModelName);
		if (!model) {
		  Rprintf("model creation failed");
		}


		/* PDF calculation */
		model->pdf(Rllsum, Rlikl, Rliku, Rllikl, Rlliku, rtl, rtu, phi);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,likl);
		SET_VECTOR_ELT(out,1,liku);
		SET_VECTOR_ELT(out,2,llikl);
		SET_VECTOR_ELT(out,3,lliku);
		SET_VECTOR_ELT(out,4,llsum);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("likl"));
		SET_STRING_ELT(names,1,Rf_mkChar("liku"));
		SET_STRING_ELT(names,2,Rf_mkChar("loglikl"));
		SET_STRING_ELT(names,3,Rf_mkChar("logliku"));
		SET_STRING_ELT(names,4,Rf_mkChar("sum_log_pdf"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		R_Free(phi);


		return(out);
	}

}



extern "C" {

  SEXP CDF(SEXP re, SEXP in, SEXP re_l, SEXP re_u, SEXP ch) {

    /* define input variables */

    ModelName = R_CHAR(STRING_ELT(ch, 0));

    N_deps = INTEGER(in)[0];
    N_rtl = INTEGER(in)[1];
    N_rtu = INTEGER(in)[2];
    dt_scale = REAL(re)[0];
    rt_max = REAL(re)[1];

    std::vector<double> rtl(N_rtl);
    std::vector<double> rtu(N_rtu);
    for (int i = 0; i < N_rtl; i++) {
      rtl[i] = REAL(re_l)[i];
    }
    for (int i = 0; i < N_rtu; i++) {
      rtu[i] = REAL(re_u)[i];
    }
    N_phi = INTEGER(in)[3];
    double *phi = (double*)R_Calloc(N_phi, double);
    for(int i=0; i<N_phi; i++) {
      phi[i] = REAL(re)[i+2];
    }

    /* declare R objects for output */
    int outCnt = 0, prtCnt = 0;
    SEXP CDFlow = PROTECT(Rf_allocVector(REALSXP, N_rtl));
    outCnt++;
    SEXP CDFupp = PROTECT(Rf_allocVector(REALSXP, N_rtu));
    outCnt++;
    SEXP logCDFlow = PROTECT(Rf_allocVector(REALSXP, N_rtl));
    outCnt++;
    SEXP logCDFupp = PROTECT(Rf_allocVector(REALSXP, N_rtu));
    outCnt++;
    SEXP sumlogCDF = PROTECT(Rf_allocVector(REALSXP, 1));
    outCnt++;
    SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
    prtCnt = outCnt + 1;


    /* declare C++ pointers for R objects */
    double *RCDFlow = REAL(CDFlow);
    double *RCDFupp = REAL(CDFupp);
    double *RlogCDFlow = REAL(logCDFlow);
    double *RlogCDFupp = REAL(logCDFupp);
    double *RsumlogCDF = REAL(sumlogCDF);


    /* model creation */
    auto model = createModel(ModelName);
    if (!model) {
      Rprintf("model creation failed");
    }


    /* CDF calculation */
    model->cdf(RsumlogCDF, RCDFlow, RCDFupp, RlogCDFlow, RlogCDFupp, rtl, rtu, phi);


    /* set elements of list out */
    SET_VECTOR_ELT(out,0,CDFlow);
    SET_VECTOR_ELT(out,1,CDFupp);
    SET_VECTOR_ELT(out,2,logCDFlow);
    SET_VECTOR_ELT(out,3,logCDFupp);
    SET_VECTOR_ELT(out,4,sumlogCDF);


    /* make name vector and set element names */
    SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
    prtCnt++;
    SET_STRING_ELT(names,0,Rf_mkChar("CDFlow"));
    SET_STRING_ELT(names,1,Rf_mkChar("CDFupp"));
    SET_STRING_ELT(names,2,Rf_mkChar("logCDFlow"));
    SET_STRING_ELT(names,3,Rf_mkChar("logCDFupp"));
    SET_STRING_ELT(names,4,Rf_mkChar("sum_log_cdf"));

    Rf_setAttrib(out,R_NamesSymbol,names);


    /* Unprotect the out and names objects */
    UNPROTECT(prtCnt);

    R_Free(phi);


    return(out);
  }

}



extern "C" {

	SEXP SIM(SEXP re, SEXP in, SEXP ch) {


		/* define input variables */

		ModelName = R_CHAR(STRING_ELT(ch, 0));

		N = INTEGER(in)[0];
		N_phi = INTEGER(in)[1];

		dt_ = REAL(re)[0];
		double *phi = (double*)R_Calloc(N_phi, double);
		for(int i=0; i<N_phi; i++) {
		  phi[i] = REAL(re)[i+1];
		}


		/* declare R objects for output */
		int outCnt = 0, prtCnt = 0;
		SEXP rt = PROTECT(Rf_allocVector(REALSXP, N));
		outCnt++;
		SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
		prtCnt = outCnt + 1;


		/* declare C++ pointers for R objects */
		double *Rrt = REAL(rt);


		/* model creation */
		auto model = createModel(ModelName);
		if (!model) {
		  Rprintf("model creation failed");
		}


		/* sampling function */
		model->rand(Rrt, phi);


		/* set elements of list out */
		SET_VECTOR_ELT(out,0,rt);


		/* make name vector and set element names */
		SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
		prtCnt++;
		SET_STRING_ELT(names,0,Rf_mkChar("rt"));

		Rf_setAttrib(out,R_NamesSymbol,names);


		/* Unprotect the out and names objects */
		UNPROTECT(prtCnt);

		R_Free(phi);


		return(out);
	}

}




extern "C" {

  SEXP grid_pdf(SEXP re, SEXP in, SEXP ch) {

    /* define input variables */
    ModelName = R_CHAR(STRING_ELT(ch, 0));

    N_deps = INTEGER(in)[0];
    N_phi = INTEGER(in)[1];

    dt_scale = REAL(re)[0];
    rt_max = REAL(re)[1];

    double *phi = (double*)R_Calloc(N_phi, double);
    for(int i=0; i<N_phi; i++) {
      phi[i] = REAL(re)[i+2];
    }

    /* declare R objects for output */
    int outCnt = 0, prtCnt = 0;
    SEXP rt = PROTECT(Rf_allocVector(REALSXP, N_dt));
    outCnt++;
    SEXP pdf_u = PROTECT(Rf_allocVector(REALSXP, N_dt));
    outCnt++;
    SEXP pdf_l = PROTECT(Rf_allocVector(REALSXP, N_dt));
    outCnt++;
    SEXP out = PROTECT(Rf_allocVector(VECSXP, outCnt));
    prtCnt = outCnt + 1;


    /* declare C++ pointers for R objects */
    double *Rrt = REAL(rt);
    double *Rpdf_u = REAL(pdf_u);
    double *Rpdf_l = REAL(pdf_l);

    /* model creation */
    auto model = createModel(ModelName);
    if (!model) {
      Rprintf("model creation failed");
    }


    /* main likelihood function */
    model->grid_pdf(Rrt, Rpdf_u, Rpdf_l, phi);

    /* set elements of list out */
    SET_VECTOR_ELT(out,0,rt);
    SET_VECTOR_ELT(out,1,pdf_u);
    SET_VECTOR_ELT(out,2,pdf_l);


    /* make name vector and set element names */
    SEXP names = PROTECT(Rf_allocVector(STRSXP, outCnt));
    prtCnt++;
    SET_STRING_ELT(names,0,Rf_mkChar("rt"));
    SET_STRING_ELT(names,1,Rf_mkChar("pdf_u"));
    SET_STRING_ELT(names,2,Rf_mkChar("pdf_l"));
    Rf_setAttrib(out,R_NamesSymbol,names);


    /* Unprotect the out and names objects */
    UNPROTECT(prtCnt);

    R_Free(phi);


    return(out);
  }

}
