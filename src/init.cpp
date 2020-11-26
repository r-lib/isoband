#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// Defined in clip-lines.cpp
extern "C" SEXP clip_lines_impl(SEXP x, SEXP y, SEXP id, SEXP _p_mid_x, SEXP _p_mid_y,
                     SEXP _width, SEXP _height, SEXP _theta, SEXP _asp);
// Defined in isobands.cpp
extern "C" SEXP isobands_impl(SEXP x, SEXP y, SEXP z, SEXP value_low, SEXP value_high);
extern "C" SEXP isolines_impl(SEXP x, SEXP y, SEXP z, SEXP value);
// Defined in separate-polygons.cpp
extern "C" SEXP separate_polygons(SEXP x, SEXP y, SEXP id);

// From testthat
extern "C" SEXP run_testthat_tests(SEXP use_xml_sxp);

static const R_CallMethodDef CallEntries[] = {
  {"clip_lines_impl_c", (DL_FUNC) &clip_lines_impl, 9},
  {"isobands_impl_c", (DL_FUNC) &isobands_impl, 5},
  {"isolines_impl_c", (DL_FUNC) &isolines_impl, 4},
  {"separate_polygons_c", (DL_FUNC) &separate_polygons, 3},
  {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 1},
  {NULL, NULL, 0}
};

extern "C" void R_init_isoband(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
