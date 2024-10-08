// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Update_logit
List Update_logit(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v, int max_t, int c, int r);
RcppExport SEXP _DiscreteKL_Update_logit(SEXP tSEXP, SEXP XSEXP, SEXP indSEXP, SEXP beta_tSEXP, SEXP beta_vSEXP, SEXP max_tSEXP, SEXP cSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type t(tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type ind(indSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v(beta_vSEXP);
    Rcpp::traits::input_parameter< int >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(Update_logit(t, X, ind, beta_t, beta_v, max_t, c, r));
    return rcpp_result_gen;
END_RCPP
}
// UpdateKL_logit
List UpdateKL_logit(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd delta, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v, int max_t, int c, int r, Eigen::MatrixXd beta_t_tilde, Eigen::MatrixXd beta_v_tilde, double eta);
RcppExport SEXP _DiscreteKL_UpdateKL_logit(SEXP tSEXP, SEXP XSEXP, SEXP deltaSEXP, SEXP beta_tSEXP, SEXP beta_vSEXP, SEXP max_tSEXP, SEXP cSEXP, SEXP rSEXP, SEXP beta_t_tildeSEXP, SEXP beta_v_tildeSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type t(tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v(beta_vSEXP);
    Rcpp::traits::input_parameter< int >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t_tilde(beta_t_tildeSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v_tilde(beta_v_tildeSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateKL_logit(t, X, delta, beta_t, beta_v, max_t, c, r, beta_t_tilde, beta_v_tilde, eta));
    return rcpp_result_gen;
END_RCPP
}
// NR_logit
List NR_logit(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v, double tol, int max_iter);
RcppExport SEXP _DiscreteKL_NR_logit(SEXP tSEXP, SEXP XSEXP, SEXP indSEXP, SEXP beta_tSEXP, SEXP beta_vSEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type t(tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type ind(indSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v(beta_vSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(NR_logit(t, X, ind, beta_t, beta_v, tol, max_iter));
    return rcpp_result_gen;
END_RCPP
}
// NRKL_logit
List NRKL_logit(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v, double tol, int max_iter, Eigen::MatrixXd beta_t_tilde, Eigen::MatrixXd beta_v_tilde, double eta);
RcppExport SEXP _DiscreteKL_NRKL_logit(SEXP tSEXP, SEXP XSEXP, SEXP indSEXP, SEXP beta_tSEXP, SEXP beta_vSEXP, SEXP tolSEXP, SEXP max_iterSEXP, SEXP beta_t_tildeSEXP, SEXP beta_v_tildeSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type t(tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type ind(indSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v(beta_vSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t_tilde(beta_t_tildeSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v_tilde(beta_v_tildeSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(NRKL_logit(t, X, ind, beta_t, beta_v, tol, max_iter, beta_t_tilde, beta_v_tilde, eta));
    return rcpp_result_gen;
END_RCPP
}
// UpdateKL2_logit
List UpdateKL2_logit(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd delta, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v, int max_t, int c, int r, Eigen::MatrixXd beta_t_tilde1, Eigen::MatrixXd beta_v_tilde1, double eta1, Eigen::MatrixXd beta_t_tilde2, Eigen::MatrixXd beta_v_tilde2, double eta2);
RcppExport SEXP _DiscreteKL_UpdateKL2_logit(SEXP tSEXP, SEXP XSEXP, SEXP deltaSEXP, SEXP beta_tSEXP, SEXP beta_vSEXP, SEXP max_tSEXP, SEXP cSEXP, SEXP rSEXP, SEXP beta_t_tilde1SEXP, SEXP beta_v_tilde1SEXP, SEXP eta1SEXP, SEXP beta_t_tilde2SEXP, SEXP beta_v_tilde2SEXP, SEXP eta2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type t(tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v(beta_vSEXP);
    Rcpp::traits::input_parameter< int >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t_tilde1(beta_t_tilde1SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v_tilde1(beta_v_tilde1SEXP);
    Rcpp::traits::input_parameter< double >::type eta1(eta1SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t_tilde2(beta_t_tilde2SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v_tilde2(beta_v_tilde2SEXP);
    Rcpp::traits::input_parameter< double >::type eta2(eta2SEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateKL2_logit(t, X, delta, beta_t, beta_v, max_t, c, r, beta_t_tilde1, beta_v_tilde1, eta1, beta_t_tilde2, beta_v_tilde2, eta2));
    return rcpp_result_gen;
END_RCPP
}
// NRKL2_logit
List NRKL2_logit(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v, double tol, int max_iter, Eigen::MatrixXd beta_t_tilde1, Eigen::MatrixXd beta_v_tilde1, double eta1, Eigen::MatrixXd beta_t_tilde2, Eigen::MatrixXd beta_v_tilde2, double eta2);
RcppExport SEXP _DiscreteKL_NRKL2_logit(SEXP tSEXP, SEXP XSEXP, SEXP indSEXP, SEXP beta_tSEXP, SEXP beta_vSEXP, SEXP tolSEXP, SEXP max_iterSEXP, SEXP beta_t_tilde1SEXP, SEXP beta_v_tilde1SEXP, SEXP eta1SEXP, SEXP beta_t_tilde2SEXP, SEXP beta_v_tilde2SEXP, SEXP eta2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type t(tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type ind(indSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v(beta_vSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t_tilde1(beta_t_tilde1SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v_tilde1(beta_v_tilde1SEXP);
    Rcpp::traits::input_parameter< double >::type eta1(eta1SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t_tilde2(beta_t_tilde2SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v_tilde2(beta_v_tilde2SEXP);
    Rcpp::traits::input_parameter< double >::type eta2(eta2SEXP);
    rcpp_result_gen = Rcpp::wrap(NRKL2_logit(t, X, ind, beta_t, beta_v, tol, max_iter, beta_t_tilde1, beta_v_tilde1, eta1, beta_t_tilde2, beta_v_tilde2, eta2));
    return rcpp_result_gen;
END_RCPP
}
// Update_cloglog
List Update_cloglog(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v, int max_t, int c, int r, double epsilon);
RcppExport SEXP _DiscreteKL_Update_cloglog(SEXP tSEXP, SEXP XSEXP, SEXP indSEXP, SEXP beta_tSEXP, SEXP beta_vSEXP, SEXP max_tSEXP, SEXP cSEXP, SEXP rSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type t(tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type ind(indSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v(beta_vSEXP);
    Rcpp::traits::input_parameter< int >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(Update_cloglog(t, X, ind, beta_t, beta_v, max_t, c, r, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// UpdateKL_cloglog
List UpdateKL_cloglog(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd delta, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v, int max_t, int c, int r, Eigen::MatrixXd beta_t_tilde, Eigen::MatrixXd beta_v_tilde, double eta, double epsilon);
RcppExport SEXP _DiscreteKL_UpdateKL_cloglog(SEXP tSEXP, SEXP XSEXP, SEXP deltaSEXP, SEXP beta_tSEXP, SEXP beta_vSEXP, SEXP max_tSEXP, SEXP cSEXP, SEXP rSEXP, SEXP beta_t_tildeSEXP, SEXP beta_v_tildeSEXP, SEXP etaSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type t(tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v(beta_vSEXP);
    Rcpp::traits::input_parameter< int >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t_tilde(beta_t_tildeSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v_tilde(beta_v_tildeSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateKL_cloglog(t, X, delta, beta_t, beta_v, max_t, c, r, beta_t_tilde, beta_v_tilde, eta, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// NR_cloglog
List NR_cloglog(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v, double tol, int max_iter, double epsilon);
RcppExport SEXP _DiscreteKL_NR_cloglog(SEXP tSEXP, SEXP XSEXP, SEXP indSEXP, SEXP beta_tSEXP, SEXP beta_vSEXP, SEXP tolSEXP, SEXP max_iterSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type t(tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type ind(indSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v(beta_vSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(NR_cloglog(t, X, ind, beta_t, beta_v, tol, max_iter, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// NRKL_cloglog
List NRKL_cloglog(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v, double tol, int max_iter, Eigen::MatrixXd beta_t_tilde, Eigen::MatrixXd beta_v_tilde, double eta, double epsilon);
RcppExport SEXP _DiscreteKL_NRKL_cloglog(SEXP tSEXP, SEXP XSEXP, SEXP indSEXP, SEXP beta_tSEXP, SEXP beta_vSEXP, SEXP tolSEXP, SEXP max_iterSEXP, SEXP beta_t_tildeSEXP, SEXP beta_v_tildeSEXP, SEXP etaSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type t(tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type ind(indSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t(beta_tSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v(beta_vSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_t_tilde(beta_t_tildeSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type beta_v_tilde(beta_v_tildeSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(NRKL_cloglog(t, X, ind, beta_t, beta_v, tol, max_iter, beta_t_tilde, beta_v_tilde, eta, epsilon));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DiscreteKL_Update_logit", (DL_FUNC) &_DiscreteKL_Update_logit, 8},
    {"_DiscreteKL_UpdateKL_logit", (DL_FUNC) &_DiscreteKL_UpdateKL_logit, 11},
    {"_DiscreteKL_NR_logit", (DL_FUNC) &_DiscreteKL_NR_logit, 7},
    {"_DiscreteKL_NRKL_logit", (DL_FUNC) &_DiscreteKL_NRKL_logit, 10},
    {"_DiscreteKL_UpdateKL2_logit", (DL_FUNC) &_DiscreteKL_UpdateKL2_logit, 14},
    {"_DiscreteKL_NRKL2_logit", (DL_FUNC) &_DiscreteKL_NRKL2_logit, 13},
    {"_DiscreteKL_Update_cloglog", (DL_FUNC) &_DiscreteKL_Update_cloglog, 9},
    {"_DiscreteKL_UpdateKL_cloglog", (DL_FUNC) &_DiscreteKL_UpdateKL_cloglog, 12},
    {"_DiscreteKL_NR_cloglog", (DL_FUNC) &_DiscreteKL_NR_cloglog, 8},
    {"_DiscreteKL_NRKL_cloglog", (DL_FUNC) &_DiscreteKL_NRKL_cloglog, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_DiscreteKL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
