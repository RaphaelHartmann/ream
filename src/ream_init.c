
#include <stdlib.h> // for NULL
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP grid_pdf(SEXP, SEXP, SEXP);
extern SEXP PDF(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP CDF(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP SIM(SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    {"grid_pdf", (DL_FUNC) &grid_pdf, 3},
    {"PDF", (DL_FUNC) &PDF, 5},
    {"CDF", (DL_FUNC) &CDF, 5},
    {"SIM", (DL_FUNC) &SIM, 3},

    {NULL, NULL, 0}
};

void R_init_ream(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
