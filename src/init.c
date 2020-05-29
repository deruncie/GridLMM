#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _GridLMM_build_downdate_Xs(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GridLMM_Calculate_qt_LASSO(SEXP, SEXP, SEXP);
extern SEXP _GridLMM_chol_c(SEXP);
extern SEXP _GridLMM_chol_dropRows(SEXP, SEXP, SEXP);
extern SEXP _GridLMM_chol_update(SEXP, SEXP, SEXP);
extern SEXP _GridLMM_chol_update_L(SEXP, SEXP, SEXP);
extern SEXP _GridLMM_crossprod_cholR(SEXP, SEXP);
extern SEXP _GridLMM_F_hats(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GridLMM_GridLMM_SS_downdate_matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GridLMM_GridLMM_SS_matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GridLMM_log_det_of_XtX(SEXP, SEXP, SEXP);
extern SEXP _GridLMM_premultiply_list_of_matrices(SEXP, SEXP);
extern SEXP _GridLMM_svd_c(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_GridLMM_build_downdate_Xs",            (DL_FUNC) &_GridLMM_build_downdate_Xs,            4},
    {"_GridLMM_Calculate_qt_LASSO",           (DL_FUNC) &_GridLMM_Calculate_qt_LASSO,           3},
    {"_GridLMM_chol_c",                       (DL_FUNC) &_GridLMM_chol_c,                       1},
    {"_GridLMM_chol_dropRows",                (DL_FUNC) &_GridLMM_chol_dropRows,                3},
    {"_GridLMM_chol_update",                  (DL_FUNC) &_GridLMM_chol_update,                  3},
    {"_GridLMM_chol_update_L",                (DL_FUNC) &_GridLMM_chol_update_L,                3},
    {"_GridLMM_crossprod_cholR",              (DL_FUNC) &_GridLMM_crossprod_cholR,              2},
    {"_GridLMM_F_hats",                       (DL_FUNC) &_GridLMM_F_hats,                       6},
    {"_GridLMM_GridLMM_SS_downdate_matrix",   (DL_FUNC) &_GridLMM_GridLMM_SS_downdate_matrix,   8},
    {"_GridLMM_GridLMM_SS_matrix",            (DL_FUNC) &_GridLMM_GridLMM_SS_matrix,            6},
    {"_GridLMM_log_det_of_XtX",               (DL_FUNC) &_GridLMM_log_det_of_XtX,               3},
    {"_GridLMM_premultiply_list_of_matrices", (DL_FUNC) &_GridLMM_premultiply_list_of_matrices, 2},
    {"_GridLMM_svd_c",                        (DL_FUNC) &_GridLMM_svd_c,                        1},
    {NULL, NULL, 0}
};

void R_init_GridLMM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
