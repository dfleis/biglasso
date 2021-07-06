// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// get_eta
RcppExport SEXP get_eta(SEXP xP, SEXP row_idx_, SEXP beta, SEXP idx_p, SEXP idx_l);
RcppExport SEXP _biglasso_get_eta(SEXP xPSEXP, SEXP row_idx_SEXP, SEXP betaSEXP, SEXP idx_pSEXP, SEXP idx_lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xP(xPSEXP);
    Rcpp::traits::input_parameter< SEXP >::type row_idx_(row_idx_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type idx_p(idx_pSEXP);
    Rcpp::traits::input_parameter< SEXP >::type idx_l(idx_lSEXP);
    rcpp_result_gen = Rcpp::wrap(get_eta(xP, row_idx_, beta, idx_p, idx_l));
    return rcpp_result_gen;
END_RCPP
}
// loglik_cox
SEXP loglik_cox(SEXP xP, SEXP offset_, SEXP row_idx_, SEXP beta, SEXP idx_p, SEXP idx_l, SEXP D_dR_sets_, SEXP d_);
RcppExport SEXP _biglasso_loglik_cox(SEXP xPSEXP, SEXP offset_SEXP, SEXP row_idx_SEXP, SEXP betaSEXP, SEXP idx_pSEXP, SEXP idx_lSEXP, SEXP D_dR_sets_SEXP, SEXP d_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xP(xPSEXP);
    Rcpp::traits::input_parameter< SEXP >::type offset_(offset_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type row_idx_(row_idx_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type idx_p(idx_pSEXP);
    Rcpp::traits::input_parameter< SEXP >::type idx_l(idx_lSEXP);
    Rcpp::traits::input_parameter< SEXP >::type D_dR_sets_(D_dR_sets_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type d_(d_SEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_cox(xP, offset_, row_idx_, beta, idx_p, idx_l, D_dR_sets_, d_));
    return rcpp_result_gen;
END_RCPP
}
