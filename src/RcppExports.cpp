// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// SimulaGGs
Rcpp::List SimulaGGs(Rcpp::IntegerVector NDatos, Rcpp::NumericVector Programa, Rcpp::CharacterVector dist_atencion, Rcpp::NumericVector param_atencion, Rcpp::CharacterVector dist_llegadas, Rcpp::NumericVector param_llegadas);
RcppExport SEXP _QueueSimulator_SimulaGGs(SEXP NDatosSEXP, SEXP ProgramaSEXP, SEXP dist_atencionSEXP, SEXP param_atencionSEXP, SEXP dist_llegadasSEXP, SEXP param_llegadasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type NDatos(NDatosSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Programa(ProgramaSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type dist_atencion(dist_atencionSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param_atencion(param_atencionSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type dist_llegadas(dist_llegadasSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type param_llegadas(param_llegadasSEXP);
    rcpp_result_gen = Rcpp::wrap(SimulaGGs(NDatos, Programa, dist_atencion, param_atencion, dist_llegadas, param_llegadas));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_QueueSimulator_SimulaGGs", (DL_FUNC) &_QueueSimulator_SimulaGGs, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_QueueSimulator(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
