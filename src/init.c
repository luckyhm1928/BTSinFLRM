#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _BTSinFLRM_inferFLRM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BTSinFLRM_orthoL2Equidense(SEXP, SEXP);
extern SEXP _BTSinFLRM_PBinner1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BTSinFLRM_PBinner2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BTSinFLRM_PBouter(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BTSinFLRM_RBinner(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BTSinFLRM_RBouter(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BTSinFLRM_WBinner(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BTSinFLRM_WBouter(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BTSinFLRM_Xsim(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BTSinFLRM_Ysim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_BTSinFLRM_inferFLRM",        (DL_FUNC) &_BTSinFLRM_inferFLRM,         8},
  {"_BTSinFLRM_orthoL2Equidense", (DL_FUNC) &_BTSinFLRM_orthoL2Equidense,  2},
  {"_BTSinFLRM_PBinner1",         (DL_FUNC) &_BTSinFLRM_PBinner1,         13},
  {"_BTSinFLRM_PBinner2",         (DL_FUNC) &_BTSinFLRM_PBinner2,         12},
  {"_BTSinFLRM_PBouter",          (DL_FUNC) &_BTSinFLRM_PBouter,          20},
  {"_BTSinFLRM_RBinner",          (DL_FUNC) &_BTSinFLRM_RBinner,          18},
  {"_BTSinFLRM_RBouter",          (DL_FUNC) &_BTSinFLRM_RBouter,          21},
  {"_BTSinFLRM_WBinner",          (DL_FUNC) &_BTSinFLRM_WBinner,          18},
  {"_BTSinFLRM_WBouter",          (DL_FUNC) &_BTSinFLRM_WBouter,          21},
  {"_BTSinFLRM_Xsim",             (DL_FUNC) &_BTSinFLRM_Xsim,              5},
  {"_BTSinFLRM_Ysim",             (DL_FUNC) &_BTSinFLRM_Ysim,              8},
  {NULL, NULL, 0}
};

void R_init_BTSinFLRM(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
