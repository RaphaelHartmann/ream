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
#include "models_7p.h"
#include "tools.h"
#include <string>
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

// for DDM with 7 parameters
int tnd_dist;
int N_dtau = 100;
double tnd_range = 3.5;
int w_dist;
int v_dist;
int N_dv = 10;
double v_range = 4.0;

  // for likelihood
// const char *OUTPUT2;
// const char *OUTPUT3;




/* Definition of the createModel function */
std::unique_ptr<Model> createModel(const char* modelName) {
  std::string modelNameStr(modelName);

  if (modelNameStr == "DMC") {
    return std::make_unique<DMC>();
  } else if (modelNameStr == "DDM") {
    return std::make_unique<DDM>();
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

	  std::string ModelNameStr(ModelName);
	  if (ModelNameStr == "DDM") {
	    tnd_dist = INTEGER(in)[4];
	    w_dist = INTEGER(in)[5];
	    v_dist = INTEGER(in)[6];
	  } else{
	    v_dist = 99;
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

    std::string ModelNameStr(ModelName);
    if (ModelNameStr == "DDM") {
      tnd_dist = INTEGER(in)[4];
      w_dist = INTEGER(in)[5];
      v_dist = INTEGER(in)[6];
    } else{
      v_dist = 99;
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

		std::string ModelNameStr(ModelName);
		if (ModelNameStr == "DDM") {
		  tnd_dist = INTEGER(in)[2];
		  w_dist = INTEGER(in)[3];
		  v_dist = INTEGER(in)[4];
		} else{
		  v_dist = 99;
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

    std::string ModelNameStr(ModelName);
    tnd_dist = INTEGER(in)[2];
    w_dist = INTEGER(in)[3];
    if (ModelNameStr == "DDM") {
      v_dist = INTEGER(in)[4];
    } else{
      v_dist = 99;
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



// ------ CUSTOM FUNCTIONS ------

// helper for preserving and releasing R objects
auto preserveReplace = [](SEXP &slot, SEXP newfun) {
  if (slot != R_NilValue)
    R_ReleaseObject(slot);  // release old function
  slot = newfun;
  R_PreserveObject(slot);   // preserve the new one
};

// -- CUSTOM FUNCTION T X --
extern "C" SEXP register_callbacks_tx(SEXP method, SEXP fnptr)
{

  if (!Rf_isString(method) || LENGTH(method) != 1)
    Rf_error("method must be a single string");

  const char* name = CHAR(STRING_ELT(method, 0));

  // copy the current callbacks so unchanged ones remain
  ModelTX_Callbacks& cb = CSTM_TX::get_callbacks();

  // R function pointer, or C++ function pointer?
  if (Rf_isFunction(fnptr)) {

    if      (strcmp(name, "drift") == 0)
      preserveReplace(cb.r_drift, fnptr);
    else if (strcmp(name, "diffusion") == 0)
      preserveReplace(cb.r_diffusion, fnptr);
    else if (strcmp(name, "upper_threshold") == 0)
      preserveReplace(cb.r_upper_threshold, fnptr);
    else if (strcmp(name, "lower_threshold") == 0)
      preserveReplace(cb.r_lower_threshold, fnptr);
    else if (strcmp(name, "non_decision") == 0)
      preserveReplace(cb.r_non_decision, fnptr);
    else if (strcmp(name, "relative_start") == 0)
      preserveReplace(cb.r_relative_start, fnptr);
    else if (strcmp(name, "contamination_strength") == 0)
      preserveReplace(cb.r_contamination_strength, fnptr);
    else if (strcmp(name, "contamination_probability") == 0)
      preserveReplace(cb.r_contamination_probability, fnptr);
    else if (strcmp(name, "modify_dt") == 0)
      preserveReplace(cb.r_modify_dt, fnptr);
    else
      Rf_error("unknown R callback name: %s", name);

  } else if (TYPEOF(fnptr) == EXTPTRSXP) {

    void* p = R_ExternalPtrAddr(fnptr);

    if      (strcmp(name, "drift") == 0)
      cb.drift = reinterpret_cast<ModelTX_Callbacks::Fn3>(p);
    else if (strcmp(name, "diffusion") == 0)
      cb.diffusion = reinterpret_cast<ModelTX_Callbacks::Fn3>(p);
    else if (strcmp(name, "upper_threshold") == 0)
      cb.upper_threshold = reinterpret_cast<ModelTX_Callbacks::Fn2>(p);
    else if (strcmp(name, "lower_threshold") == 0)
      cb.lower_threshold = reinterpret_cast<ModelTX_Callbacks::Fn2>(p);
    else if (strcmp(name, "non_decision") == 0)
      cb.non_decision = reinterpret_cast<ModelTX_Callbacks::Fn1>(p);
    else if (strcmp(name, "relative_start") == 0)
      cb.relative_start = reinterpret_cast<ModelTX_Callbacks::Fn1>(p);
    else if (strcmp(name, "contamination_strength") == 0)
      cb.contamination_strength = reinterpret_cast<ModelTX_Callbacks::Fn1>(p);
    else if (strcmp(name, "contamination_probability") == 0)
      cb.contamination_probability = reinterpret_cast<ModelTX_Callbacks::Fn2>(p);
    else if (strcmp(name, "modify_dt") == 0)
      cb.modify_dt = reinterpret_cast<ModelTX_Callbacks::Fn2>(p);
    else
      Rf_error("unknown method name: %s", name);

  } else if (fnptr == R_NilValue) {
    // explicit remove
    if      (strcmp(name, "drift") == 0) cb.r_drift = R_NilValue;
    else if (strcmp(name, "diffusion") == 0) cb.r_diffusion = R_NilValue;
    else if (strcmp(name, "upper_threshold") == 0) cb.r_upper_threshold = R_NilValue;
    else if (strcmp(name, "lower_threshold") == 0) cb.r_lower_threshold = R_NilValue;
    else if (strcmp(name, "non_decision") == 0) cb.r_non_decision = R_NilValue;
    else if (strcmp(name, "relative_start") == 0) cb.r_relative_start = R_NilValue;
    else if (strcmp(name, "contamination_strength") == 0) cb.r_contamination_strength = R_NilValue;
    else if (strcmp(name, "contamination_probability") == 0) cb.r_contamination_probability = R_NilValue;
    else if (strcmp(name, "modify_dt") == 0) cb.r_modify_dt = R_NilValue;
    CSTM_TX::set_callbacks(cb);
    return R_NilValue;
  } else {
    Rf_error("fnptr must be either a function or an external pointer");
  }

  CSTM_TX::set_callbacks(cb);
  return R_NilValue;

}

extern "C" SEXP unregister_callbacks_tx()
{
  // Retrieve the current callback set
  ModelTX_Callbacks old = CSTM_TX::get_callbacks();

  // Release every preserved R function if it exists
  if (old.r_drift != R_NilValue)
    R_ReleaseObject(old.r_drift);
  if (old.r_diffusion != R_NilValue)
    R_ReleaseObject(old.r_diffusion);
  if (old.r_upper_threshold != R_NilValue)
    R_ReleaseObject(old.r_upper_threshold);
  if (old.r_lower_threshold != R_NilValue)
    R_ReleaseObject(old.r_lower_threshold);
  if (old.r_non_decision != R_NilValue)
    R_ReleaseObject(old.r_non_decision);
  if (old.r_relative_start != R_NilValue)
    R_ReleaseObject(old.r_relative_start);
  if (old.r_contamination_strength != R_NilValue)
    R_ReleaseObject(old.r_contamination_strength);
  if (old.r_contamination_probability != R_NilValue)
    R_ReleaseObject(old.r_contamination_probability);
  if (old.r_modify_dt != R_NilValue)
    R_ReleaseObject(old.r_modify_dt);

  // Replace current callbacks with an empty (all–null) set
  ModelTX_Callbacks empty;
  CSTM_TX::set_callbacks(empty);

  return R_NilValue;
}



// -- CUSTOM FUNCTION T --
extern "C" SEXP register_callbacks_t(SEXP method, SEXP fnptr)
{

  if (!Rf_isString(method) || LENGTH(method) != 1)
    Rf_error("method must be a single string");

  const char* name = CHAR(STRING_ELT(method, 0));

  // copy the current callbacks so unchanged ones remain
  ModelT_Callbacks& cb = CSTM_T::get_callbacks();

  // R function pointer, or C++ function pointer?
  if (Rf_isFunction(fnptr)) {

    if      (strcmp(name, "drift") == 0)
      preserveReplace(cb.r_drift, fnptr);
    else if (strcmp(name, "diffusion") == 0)
      preserveReplace(cb.r_diffusion, fnptr);
    else if (strcmp(name, "upper_threshold") == 0)
      preserveReplace(cb.r_upper_threshold, fnptr);
    else if (strcmp(name, "lower_threshold") == 0)
      preserveReplace(cb.r_lower_threshold, fnptr);
    else if (strcmp(name, "non_decision") == 0)
      preserveReplace(cb.r_non_decision, fnptr);
    else if (strcmp(name, "relative_start") == 0)
      preserveReplace(cb.r_relative_start, fnptr);
    else if (strcmp(name, "contamination_strength") == 0)
      preserveReplace(cb.r_contamination_strength, fnptr);
    else if (strcmp(name, "contamination_probability") == 0)
      preserveReplace(cb.r_contamination_probability, fnptr);
    else if (strcmp(name, "modify_dt") == 0)
      preserveReplace(cb.r_modify_dt, fnptr);
    else
      Rf_error("unknown R callback name: %s", name);

  } else if (TYPEOF(fnptr) == EXTPTRSXP) {

    void* p = R_ExternalPtrAddr(fnptr);

    if      (strcmp(name, "drift") == 0)
      cb.drift = reinterpret_cast<ModelT_Callbacks::Fn2>(p);
    else if (strcmp(name, "diffusion") == 0)
      cb.diffusion = reinterpret_cast<ModelT_Callbacks::Fn3>(p);
    else if (strcmp(name, "upper_threshold") == 0)
      cb.upper_threshold = reinterpret_cast<ModelT_Callbacks::Fn2>(p);
    else if (strcmp(name, "lower_threshold") == 0)
      cb.lower_threshold = reinterpret_cast<ModelT_Callbacks::Fn2>(p);
    else if (strcmp(name, "non_decision") == 0)
      cb.non_decision = reinterpret_cast<ModelT_Callbacks::Fn1>(p);
    else if (strcmp(name, "relative_start") == 0)
      cb.relative_start = reinterpret_cast<ModelT_Callbacks::Fn1>(p);
    else if (strcmp(name, "contamination_strength") == 0)
      cb.contamination_strength = reinterpret_cast<ModelT_Callbacks::Fn1>(p);
    else if (strcmp(name, "contamination_probability") == 0)
      cb.contamination_probability = reinterpret_cast<ModelT_Callbacks::Fn2>(p);
    else if (strcmp(name, "modify_dt") == 0)
      cb.modify_dt = reinterpret_cast<ModelT_Callbacks::Fn2>(p);
    else
      Rf_error("unknown method name: %s", name);

  } else if (fnptr == R_NilValue) {
    // explicit remove
    if      (strcmp(name, "drift") == 0) cb.r_drift = R_NilValue;
    else if (strcmp(name, "diffusion") == 0) cb.r_diffusion = R_NilValue;
    else if (strcmp(name, "upper_threshold") == 0) cb.r_upper_threshold = R_NilValue;
    else if (strcmp(name, "lower_threshold") == 0) cb.r_lower_threshold = R_NilValue;
    else if (strcmp(name, "non_decision") == 0) cb.r_non_decision = R_NilValue;
    else if (strcmp(name, "relative_start") == 0) cb.r_relative_start = R_NilValue;
    else if (strcmp(name, "contamination_strength") == 0) cb.r_contamination_strength = R_NilValue;
    else if (strcmp(name, "contamination_probability") == 0) cb.r_contamination_probability = R_NilValue;
    else if (strcmp(name, "modify_dt") == 0) cb.r_modify_dt = R_NilValue;
    CSTM_T::set_callbacks(cb);
    return R_NilValue;
  } else {
    Rf_error("fnptr must be either a function or an external pointer");
  }

  CSTM_T::set_callbacks(cb);
  return R_NilValue;

}

extern "C" SEXP unregister_callbacks_t()
{
  // Retrieve the current callback set
  ModelT_Callbacks old = CSTM_T::get_callbacks();

  // Release every preserved R function if it exists
  if (old.r_drift != R_NilValue)
    R_ReleaseObject(old.r_drift);
  if (old.r_diffusion != R_NilValue)
    R_ReleaseObject(old.r_diffusion);
  if (old.r_upper_threshold != R_NilValue)
    R_ReleaseObject(old.r_upper_threshold);
  if (old.r_lower_threshold != R_NilValue)
    R_ReleaseObject(old.r_lower_threshold);
  if (old.r_non_decision != R_NilValue)
    R_ReleaseObject(old.r_non_decision);
  if (old.r_relative_start != R_NilValue)
    R_ReleaseObject(old.r_relative_start);
  if (old.r_contamination_strength != R_NilValue)
    R_ReleaseObject(old.r_contamination_strength);
  if (old.r_contamination_probability != R_NilValue)
    R_ReleaseObject(old.r_contamination_probability);
  if (old.r_modify_dt != R_NilValue)
    R_ReleaseObject(old.r_modify_dt);

  // Replace current callbacks with an empty (all–null) set
  ModelT_Callbacks empty;
  CSTM_T::set_callbacks(empty);

  return R_NilValue;
}



// -- CUSTOM FUNCTION T W --
extern "C" SEXP register_callbacks_tw(SEXP method, SEXP fnptr)
{

  if (!Rf_isString(method) || LENGTH(method) != 1)
    Rf_error("method must be a single string");

  const char* name = CHAR(STRING_ELT(method, 0));

  // copy the current callbacks so unchanged ones remain
  ModelTW_Callbacks& cb = CSTM_TW::get_callbacks();

  // R function pointer, or C++ function pointer?
  if (Rf_isFunction(fnptr)) {

    if      (strcmp(name, "drift") == 0)
      preserveReplace(cb.r_drift, fnptr);
    else if (strcmp(name, "diffusion") == 0)
      preserveReplace(cb.r_diffusion, fnptr);
    else if (strcmp(name, "upper_threshold") == 0)
      preserveReplace(cb.r_upper_threshold, fnptr);
    else if (strcmp(name, "lower_threshold") == 0)
      preserveReplace(cb.r_lower_threshold, fnptr);
    else if (strcmp(name, "non_decision") == 0)
      preserveReplace(cb.r_non_decision, fnptr);
    else if (strcmp(name, "relative_start") == 0)
      preserveReplace(cb.r_relative_start, fnptr);
    else if (strcmp(name, "contamination_strength") == 0)
      preserveReplace(cb.r_contamination_strength, fnptr);
    else if (strcmp(name, "contamination_probability") == 0)
      preserveReplace(cb.r_contamination_probability, fnptr);
    else if (strcmp(name, "modify_dt") == 0)
      preserveReplace(cb.r_modify_dt, fnptr);
    else
      Rf_error("unknown R callback name: %s", name);

  } else if (TYPEOF(fnptr) == EXTPTRSXP) {

    void* p = R_ExternalPtrAddr(fnptr);

    if      (strcmp(name, "drift") == 0)
      cb.drift = reinterpret_cast<ModelTW_Callbacks::Fn3>(p);
    else if (strcmp(name, "drift_ts") == 0)
      cb.drift_ts = reinterpret_cast<ModelTW_Callbacks::Fn1>(p);
    else if (strcmp(name, "diffusion") == 0)
      cb.diffusion = reinterpret_cast<ModelTW_Callbacks::Fn3>(p);
    else if (strcmp(name, "diffusion_ts") == 0)
      cb.diffusion_ts = reinterpret_cast<ModelTW_Callbacks::Fn1>(p);
    else if (strcmp(name, "upper_threshold") == 0)
      cb.upper_threshold = reinterpret_cast<ModelTW_Callbacks::Fn2>(p);
    else if (strcmp(name, "upper_threshold_ts") == 0)
      cb.upper_threshold_ts = reinterpret_cast<ModelTW_Callbacks::Fn1>(p);
    else if (strcmp(name, "lower_threshold") == 0)
      cb.lower_threshold = reinterpret_cast<ModelTW_Callbacks::Fn2>(p);
    else if (strcmp(name, "lower_threshold_ts") == 0)
      cb.lower_threshold_ts = reinterpret_cast<ModelTW_Callbacks::Fn1>(p);
    else if (strcmp(name, "non_decision") == 0)
      cb.non_decision = reinterpret_cast<ModelTW_Callbacks::Fn1>(p);
    else if (strcmp(name, "relative_start") == 0)
      cb.relative_start = reinterpret_cast<ModelTW_Callbacks::Fn1>(p);
    else if (strcmp(name, "relative_start_ts") == 0)
      cb.relative_start_ts = reinterpret_cast<ModelTW_Callbacks::Fn1>(p);
    else if (strcmp(name, "contamination_strength") == 0)
      cb.contamination_strength = reinterpret_cast<ModelTW_Callbacks::Fn1>(p);
    else if (strcmp(name, "contamination_probability") == 0)
      cb.contamination_probability = reinterpret_cast<ModelTW_Callbacks::Fn2>(p);
    else if (strcmp(name, "modify_dt") == 0)
      cb.modify_dt = reinterpret_cast<ModelTW_Callbacks::Fn2>(p);
    else if (strcmp(name, "ts_cdf") == 0)
      cb.ts_cdf = reinterpret_cast<ModelTW_Callbacks::Fn2>(p);
    else
      Rf_error("unknown method name: %s", name);

  } else if (fnptr == R_NilValue) {
    // explicit remove
    if      (strcmp(name, "drift") == 0) cb.r_drift = R_NilValue;
    else if (strcmp(name, "drift_ts") == 0) cb.r_drift_ts = R_NilValue;
    else if (strcmp(name, "diffusion") == 0) cb.r_diffusion = R_NilValue;
    else if (strcmp(name, "diffusion_ts") == 0) cb.r_diffusion_ts = R_NilValue;
    else if (strcmp(name, "upper_threshold") == 0) cb.r_upper_threshold = R_NilValue;
    else if (strcmp(name, "upper_threshold_ts") == 0) cb.r_upper_threshold_ts = R_NilValue;
    else if (strcmp(name, "lower_threshold") == 0) cb.r_lower_threshold = R_NilValue;
    else if (strcmp(name, "lower_threshold_ts") == 0) cb.r_lower_threshold_ts = R_NilValue;
    else if (strcmp(name, "non_decision") == 0) cb.r_non_decision = R_NilValue;
    else if (strcmp(name, "relative_start") == 0) cb.r_relative_start = R_NilValue;
    else if (strcmp(name, "relative_start_ts") == 0) cb.r_relative_start_ts = R_NilValue;
    else if (strcmp(name, "contamination_strength") == 0) cb.r_contamination_strength = R_NilValue;
    else if (strcmp(name, "contamination_probability") == 0) cb.r_contamination_probability = R_NilValue;
    else if (strcmp(name, "modify_dt") == 0) cb.r_modify_dt = R_NilValue;
    else if (strcmp(name, "ts_cdf") == 0) cb.r_ts_cdf = R_NilValue;
    CSTM_TW::set_callbacks(cb);
    return R_NilValue;
  } else {
    Rf_error("fnptr must be either a function or an external pointer");
  }

  CSTM_TW::set_callbacks(cb);
  return R_NilValue;

}

extern "C" SEXP unregister_callbacks_tw()
{
  // Retrieve the current callback set
  ModelTW_Callbacks old = CSTM_TW::get_callbacks();

  // Release every preserved R function if it exists
  if (old.r_drift != R_NilValue)
    R_ReleaseObject(old.r_drift);
  if (old.r_drift_ts != R_NilValue)
    R_ReleaseObject(old.r_drift_ts);
  if (old.r_diffusion != R_NilValue)
    R_ReleaseObject(old.r_diffusion);
  if (old.r_diffusion_ts != R_NilValue)
    R_ReleaseObject(old.r_diffusion_ts);
  if (old.r_upper_threshold != R_NilValue)
    R_ReleaseObject(old.r_upper_threshold);
  if (old.r_upper_threshold_ts != R_NilValue)
    R_ReleaseObject(old.r_upper_threshold_ts);
  if (old.r_lower_threshold != R_NilValue)
    R_ReleaseObject(old.r_lower_threshold);
  if (old.r_lower_threshold_ts != R_NilValue)
    R_ReleaseObject(old.r_lower_threshold_ts);
  if (old.r_non_decision != R_NilValue)
    R_ReleaseObject(old.r_non_decision);
  if (old.r_relative_start != R_NilValue)
    R_ReleaseObject(old.r_relative_start);
  if (old.r_relative_start_ts != R_NilValue)
    R_ReleaseObject(old.r_relative_start_ts);
  if (old.r_contamination_strength != R_NilValue)
    R_ReleaseObject(old.r_contamination_strength);
  if (old.r_contamination_probability != R_NilValue)
    R_ReleaseObject(old.r_contamination_probability);
  if (old.r_modify_dt != R_NilValue)
    R_ReleaseObject(old.r_modify_dt);
  if (old.r_ts_cdf != R_NilValue)
    R_ReleaseObject(old.r_ts_cdf);

  // Replace current callbacks with an empty (all–null) set
  ModelTW_Callbacks empty;
  CSTM_TW::set_callbacks(empty);

  return R_NilValue;
}
