// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cppoverlaps
Rcpp::DataFrame cppoverlaps(const Rcpp::DataFrame& df1, const Rcpp::DataFrame& df2, bool verbose, bool index_only);
RcppExport SEXP _roverlaps_cppoverlaps(SEXP df1SEXP, SEXP df2SEXP, SEXP verboseSEXP, SEXP index_onlySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df1(df1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df2(df2SEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type index_only(index_onlySEXP);
    rcpp_result_gen = Rcpp::wrap(cppoverlaps(df1, df2, verbose, index_only));
    return rcpp_result_gen;
END_RCPP
}
// cppraggeddiff
Rcpp::NumericVector cppraggeddiff(const Rcpp::NumericVector& query, const Rcpp::NumericVector& subject, bool max, int sign);
RcppExport SEXP _roverlaps_cppraggeddiff(SEXP querySEXP, SEXP subjectSEXP, SEXP maxSEXP, SEXP signSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type query(querySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type subject(subjectSEXP);
    Rcpp::traits::input_parameter< bool >::type max(maxSEXP);
    Rcpp::traits::input_parameter< int >::type sign(signSEXP);
    rcpp_result_gen = Rcpp::wrap(cppraggeddiff(query, subject, max, sign));
    return rcpp_result_gen;
END_RCPP
}
// cppraggeddiffinterval
Rcpp::NumericVector cppraggeddiffinterval(const Rcpp::NumericVector& query, const Rcpp::NumericVector& subject_start, const Rcpp::NumericVector& subject_end, bool max, int sign);
RcppExport SEXP _roverlaps_cppraggeddiffinterval(SEXP querySEXP, SEXP subject_startSEXP, SEXP subject_endSEXP, SEXP maxSEXP, SEXP signSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type query(querySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type subject_start(subject_startSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type subject_end(subject_endSEXP);
    Rcpp::traits::input_parameter< bool >::type max(maxSEXP);
    Rcpp::traits::input_parameter< int >::type sign(signSEXP);
    rcpp_result_gen = Rcpp::wrap(cppraggeddiffinterval(query, subject_start, subject_end, max, sign));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_roverlaps_cppoverlaps", (DL_FUNC) &_roverlaps_cppoverlaps, 4},
    {"_roverlaps_cppraggeddiff", (DL_FUNC) &_roverlaps_cppraggeddiff, 4},
    {"_roverlaps_cppraggeddiffinterval", (DL_FUNC) &_roverlaps_cppraggeddiffinterval, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_roverlaps(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
