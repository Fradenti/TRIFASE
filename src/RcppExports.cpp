// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// AlgoAll_cpp
Rcpp::List AlgoAll_cpp(arma::mat X, arma::mat Fm, arma::mat Gm, arma::mat L, arma::mat Linv, arma::mat Dinv, int K, int R, int niter, int swipes, double threshold, double nugget, std::string upd_row, std::string upd_col);
RcppExport SEXP _TRIFASE_AlgoAll_cpp(SEXP XSEXP, SEXP FmSEXP, SEXP GmSEXP, SEXP LSEXP, SEXP LinvSEXP, SEXP DinvSEXP, SEXP KSEXP, SEXP RSEXP, SEXP niterSEXP, SEXP swipesSEXP, SEXP thresholdSEXP, SEXP nuggetSEXP, SEXP upd_rowSEXP, SEXP upd_colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fm(FmSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gm(GmSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Linv(LinvSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dinv(DinvSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type swipes(swipesSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< std::string >::type upd_row(upd_rowSEXP);
    Rcpp::traits::input_parameter< std::string >::type upd_col(upd_colSEXP);
    rcpp_result_gen = Rcpp::wrap(AlgoAll_cpp(X, Fm, Gm, L, Linv, Dinv, K, R, niter, swipes, threshold, nugget, upd_row, upd_col));
    return rcpp_result_gen;
END_RCPP
}
// Algo4_cpp
Rcpp::List Algo4_cpp(arma::mat X, arma::mat Fm, arma::mat Gm, arma::mat L, arma::mat Linv, arma::mat Dinv, int K, int R, int niter, int swipes, double threshold, double nugget, bool stoc);
RcppExport SEXP _TRIFASE_Algo4_cpp(SEXP XSEXP, SEXP FmSEXP, SEXP GmSEXP, SEXP LSEXP, SEXP LinvSEXP, SEXP DinvSEXP, SEXP KSEXP, SEXP RSEXP, SEXP niterSEXP, SEXP swipesSEXP, SEXP thresholdSEXP, SEXP nuggetSEXP, SEXP stocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fm(FmSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gm(GmSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Linv(LinvSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dinv(DinvSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type swipes(swipesSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< bool >::type stoc(stocSEXP);
    rcpp_result_gen = Rcpp::wrap(Algo4_cpp(X, Fm, Gm, L, Linv, Dinv, K, R, niter, swipes, threshold, nugget, stoc));
    return rcpp_result_gen;
END_RCPP
}
// aa
arma::vec aa(double x);
RcppExport SEXP _TRIFASE_aa(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(aa(x));
    return rcpp_result_gen;
END_RCPP
}
// AlgoAll_cpp2
Rcpp::List AlgoAll_cpp2(arma::mat X, arma::mat Fm, arma::mat Gm, arma::mat L, arma::mat iKernel, double current_tau, arma::mat Linv, arma::mat Dinv, int K, int R, int niter, int swipes, double threshold, double nugget, bool estimate_tau, std::string upd_row, std::string upd_col);
RcppExport SEXP _TRIFASE_AlgoAll_cpp2(SEXP XSEXP, SEXP FmSEXP, SEXP GmSEXP, SEXP LSEXP, SEXP iKernelSEXP, SEXP current_tauSEXP, SEXP LinvSEXP, SEXP DinvSEXP, SEXP KSEXP, SEXP RSEXP, SEXP niterSEXP, SEXP swipesSEXP, SEXP thresholdSEXP, SEXP nuggetSEXP, SEXP estimate_tauSEXP, SEXP upd_rowSEXP, SEXP upd_colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fm(FmSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Gm(GmSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type iKernel(iKernelSEXP);
    Rcpp::traits::input_parameter< double >::type current_tau(current_tauSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Linv(LinvSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dinv(DinvSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type swipes(swipesSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< bool >::type estimate_tau(estimate_tauSEXP);
    Rcpp::traits::input_parameter< std::string >::type upd_row(upd_rowSEXP);
    Rcpp::traits::input_parameter< std::string >::type upd_col(upd_colSEXP);
    rcpp_result_gen = Rcpp::wrap(AlgoAll_cpp2(X, Fm, Gm, L, iKernel, current_tau, Linv, Dinv, K, R, niter, swipes, threshold, nugget, estimate_tau, upd_row, upd_col));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TRIFASE_AlgoAll_cpp", (DL_FUNC) &_TRIFASE_AlgoAll_cpp, 14},
    {"_TRIFASE_Algo4_cpp", (DL_FUNC) &_TRIFASE_Algo4_cpp, 13},
    {"_TRIFASE_aa", (DL_FUNC) &_TRIFASE_aa, 1},
    {"_TRIFASE_AlgoAll_cpp2", (DL_FUNC) &_TRIFASE_AlgoAll_cpp2, 17},
    {NULL, NULL, 0}
};

RcppExport void R_init_TRIFASE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
