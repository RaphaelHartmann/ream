
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
extern SEXP register_callbacks_tx(SEXP, SEXP);
extern SEXP unregister_callbacks_tx(void);
extern SEXP register_callbacks_t(SEXP, SEXP);
extern SEXP unregister_callbacks_t(void);
extern SEXP register_callbacks_tw(SEXP, SEXP);
extern SEXP unregister_callbacks_tw(void);


static const R_CallMethodDef CallEntries[] = {
    {"grid_pdf", (DL_FUNC) &grid_pdf, 3},
    {"PDF", (DL_FUNC) &PDF, 5},
    {"CDF", (DL_FUNC) &CDF, 5},
    {"SIM", (DL_FUNC) &SIM, 3},
    {"register_callbacks_tx", (DL_FUNC) &register_callbacks_tx, 2},
    {"unregister_callbacks_tx", (DL_FUNC) &unregister_callbacks_tx, 0},
    {"register_callbacks_t", (DL_FUNC) &register_callbacks_t, 2},
    {"unregister_callbacks_t", (DL_FUNC) &unregister_callbacks_t, 0},
    {"register_callbacks_tw", (DL_FUNC) &register_callbacks_tw, 2},
    {"unregister_callbacks_tw", (DL_FUNC) &unregister_callbacks_tw, 0},

    {NULL, NULL, 0}
};

void R_init_ream(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
