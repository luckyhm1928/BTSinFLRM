// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// orthoL2Equidense
arma::mat orthoL2Equidense(arma::mat& X, arma::vec& tGrid);
RcppExport SEXP _BTSinFLRM_orthoL2Equidense(SEXP XSEXP, SEXP tGridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tGrid(tGridSEXP);
    rcpp_result_gen = Rcpp::wrap(orthoL2Equidense(X, tGrid));
    return rcpp_result_gen;
END_RCPP
}
// Xsim
arma::mat Xsim(int n, arma::vec& ev, arma::mat& ef, int xi_type, double nu_xi);
RcppExport SEXP _BTSinFLRM_Xsim(SEXP nSEXP, SEXP evSEXP, SEXP efSEXP, SEXP xi_typeSEXP, SEXP nu_xiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type ev(evSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ef(efSEXP);
    Rcpp::traits::input_parameter< int >::type xi_type(xi_typeSEXP);
    Rcpp::traits::input_parameter< double >::type nu_xi(nu_xiSEXP);
    rcpp_result_gen = Rcpp::wrap(Xsim(n, ev, ef, xi_type, nu_xi));
    return rcpp_result_gen;
END_RCPP
}
// Ysim
List Ysim(arma::mat& X, arma::vec& beta, arma::vec& ev, arma::mat& ef, arma::vec& tGrid, int distn_type, int var_type, int dep_type);
RcppExport SEXP _BTSinFLRM_Ysim(SEXP XSEXP, SEXP betaSEXP, SEXP evSEXP, SEXP efSEXP, SEXP tGridSEXP, SEXP distn_typeSEXP, SEXP var_typeSEXP, SEXP dep_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type ev(evSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ef(efSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tGrid(tGridSEXP);
    Rcpp::traits::input_parameter< int >::type distn_type(distn_typeSEXP);
    Rcpp::traits::input_parameter< int >::type var_type(var_typeSEXP);
    Rcpp::traits::input_parameter< int >::type dep_type(dep_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(Ysim(X, beta, ev, ef, tGrid, distn_type, var_type, dep_type));
    return rcpp_result_gen;
END_RCPP
}
// inferFLRM
List inferFLRM(arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec, double ap);
RcppExport SEXP _BTSinFLRM_inferFLRM(SEXP XSEXP, SEXP YSEXP, SEXP X0SEXP, SEXP tGridSEXP, SEXP kn_vecSEXP, SEXP hn_vecSEXP, SEXP gn_vecSEXP, SEXP apSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tGrid(tGridSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type kn_vec(kn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type hn_vec(hn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gn_vec(gn_vecSEXP);
    Rcpp::traits::input_parameter< double >::type ap(apSEXP);
    rcpp_result_gen = Rcpp::wrap(inferFLRM(X, Y, X0, tGrid, kn_vec, hn_vec, gn_vec, ap));
    return rcpp_result_gen;
END_RCPP
}
// RBinner
List RBinner(arma::vec& idx, arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, arma::mat& betaHat_mat, arma::mat& Y0Hat_mat, arma::cube& Y0Hat_trunc_hg_cube, arma::mat& erHatCent_mat, arma::mat& YHat_mat, arma::mat& XCent, arma::mat& YTilde_mat, arma::mat& Y0Tilde_mat, arma::cube& GaInvHat_cube, arma::mat& txHat_mat_h, arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec);
RcppExport SEXP _BTSinFLRM_RBinner(SEXP idxSEXP, SEXP XSEXP, SEXP YSEXP, SEXP X0SEXP, SEXP tGridSEXP, SEXP betaHat_matSEXP, SEXP Y0Hat_matSEXP, SEXP Y0Hat_trunc_hg_cubeSEXP, SEXP erHatCent_matSEXP, SEXP YHat_matSEXP, SEXP XCentSEXP, SEXP YTilde_matSEXP, SEXP Y0Tilde_matSEXP, SEXP GaInvHat_cubeSEXP, SEXP txHat_mat_hSEXP, SEXP kn_vecSEXP, SEXP hn_vecSEXP, SEXP gn_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tGrid(tGridSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type betaHat_mat(betaHat_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y0Hat_mat(Y0Hat_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Y0Hat_trunc_hg_cube(Y0Hat_trunc_hg_cubeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type erHatCent_mat(erHatCent_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type YHat_mat(YHat_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type XCent(XCentSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type YTilde_mat(YTilde_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y0Tilde_mat(Y0Tilde_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type GaInvHat_cube(GaInvHat_cubeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type txHat_mat_h(txHat_mat_hSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type kn_vec(kn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type hn_vec(hn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gn_vec(gn_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(RBinner(idx, X, Y, X0, tGrid, betaHat_mat, Y0Hat_mat, Y0Hat_trunc_hg_cube, erHatCent_mat, YHat_mat, XCent, YTilde_mat, Y0Tilde_mat, GaInvHat_cube, txHat_mat_h, kn_vec, hn_vec, gn_vec));
    return rcpp_result_gen;
END_RCPP
}
// RBouter
List RBouter(arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, arma::mat& betaHat_mat, arma::mat& Y0Hat_mat, arma::cube& Y0Hat_trunc_hg_cube, arma::mat& erHatCent_mat, arma::mat& YHat_mat, arma::mat& XCent, arma::mat& YTilde_mat, arma::mat& Y0Tilde_mat, arma::cube& GaInvHat_cube, arma::mat& txHat_mat_h, arma::mat& stat_nonstd_mat, arma::mat& stat_Tstd_mat, arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec, int M_bts, double ap);
RcppExport SEXP _BTSinFLRM_RBouter(SEXP XSEXP, SEXP YSEXP, SEXP X0SEXP, SEXP tGridSEXP, SEXP betaHat_matSEXP, SEXP Y0Hat_matSEXP, SEXP Y0Hat_trunc_hg_cubeSEXP, SEXP erHatCent_matSEXP, SEXP YHat_matSEXP, SEXP XCentSEXP, SEXP YTilde_matSEXP, SEXP Y0Tilde_matSEXP, SEXP GaInvHat_cubeSEXP, SEXP txHat_mat_hSEXP, SEXP stat_nonstd_matSEXP, SEXP stat_Tstd_matSEXP, SEXP kn_vecSEXP, SEXP hn_vecSEXP, SEXP gn_vecSEXP, SEXP M_btsSEXP, SEXP apSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tGrid(tGridSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type betaHat_mat(betaHat_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y0Hat_mat(Y0Hat_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Y0Hat_trunc_hg_cube(Y0Hat_trunc_hg_cubeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type erHatCent_mat(erHatCent_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type YHat_mat(YHat_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type XCent(XCentSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type YTilde_mat(YTilde_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y0Tilde_mat(Y0Tilde_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type GaInvHat_cube(GaInvHat_cubeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type txHat_mat_h(txHat_mat_hSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type stat_nonstd_mat(stat_nonstd_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type stat_Tstd_mat(stat_Tstd_matSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type kn_vec(kn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type hn_vec(hn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gn_vec(gn_vecSEXP);
    Rcpp::traits::input_parameter< int >::type M_bts(M_btsSEXP);
    Rcpp::traits::input_parameter< double >::type ap(apSEXP);
    rcpp_result_gen = Rcpp::wrap(RBouter(X, Y, X0, tGrid, betaHat_mat, Y0Hat_mat, Y0Hat_trunc_hg_cube, erHatCent_mat, YHat_mat, XCent, YTilde_mat, Y0Tilde_mat, GaInvHat_cube, txHat_mat_h, stat_nonstd_mat, stat_Tstd_mat, kn_vec, hn_vec, gn_vec, M_bts, ap));
    return rcpp_result_gen;
END_RCPP
}
// PBinner1
List PBinner1(arma::vec& idx, arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, arma::mat& betaHat_mat, arma::mat& UHat_mat, arma::cube& sxHat_cube_kh, arma::mat& Y0Hat_mat, arma::cube& Y0Hat_trunc_hg_cube, arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec);
RcppExport SEXP _BTSinFLRM_PBinner1(SEXP idxSEXP, SEXP XSEXP, SEXP YSEXP, SEXP X0SEXP, SEXP tGridSEXP, SEXP betaHat_matSEXP, SEXP UHat_matSEXP, SEXP sxHat_cube_khSEXP, SEXP Y0Hat_matSEXP, SEXP Y0Hat_trunc_hg_cubeSEXP, SEXP kn_vecSEXP, SEXP hn_vecSEXP, SEXP gn_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tGrid(tGridSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type betaHat_mat(betaHat_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type UHat_mat(UHat_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type sxHat_cube_kh(sxHat_cube_khSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y0Hat_mat(Y0Hat_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Y0Hat_trunc_hg_cube(Y0Hat_trunc_hg_cubeSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type kn_vec(kn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type hn_vec(hn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gn_vec(gn_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(PBinner1(idx, X, Y, X0, tGrid, betaHat_mat, UHat_mat, sxHat_cube_kh, Y0Hat_mat, Y0Hat_trunc_hg_cube, kn_vec, hn_vec, gn_vec));
    return rcpp_result_gen;
END_RCPP
}
// PBinner2
List PBinner2(arma::vec& idxTilde, arma::mat& X, arma::mat& YTilde_mat, arma::mat& X0, arma::vec& tGrid, arma::mat& betaTilde_mat, arma::mat& UTilde_mat, arma::cube& sxHat_cube_kh, arma::mat& Y0Tilde_mat, arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec);
RcppExport SEXP _BTSinFLRM_PBinner2(SEXP idxTildeSEXP, SEXP XSEXP, SEXP YTilde_matSEXP, SEXP X0SEXP, SEXP tGridSEXP, SEXP betaTilde_matSEXP, SEXP UTilde_matSEXP, SEXP sxHat_cube_khSEXP, SEXP Y0Tilde_matSEXP, SEXP kn_vecSEXP, SEXP hn_vecSEXP, SEXP gn_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type idxTilde(idxTildeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type YTilde_mat(YTilde_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tGrid(tGridSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type betaTilde_mat(betaTilde_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type UTilde_mat(UTilde_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type sxHat_cube_kh(sxHat_cube_khSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y0Tilde_mat(Y0Tilde_matSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type kn_vec(kn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type hn_vec(hn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gn_vec(gn_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(PBinner2(idxTilde, X, YTilde_mat, X0, tGrid, betaTilde_mat, UTilde_mat, sxHat_cube_kh, Y0Tilde_mat, kn_vec, hn_vec, gn_vec));
    return rcpp_result_gen;
END_RCPP
}
// PBouter
List PBouter(arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, arma::mat& betaHat_mat, arma::mat& UHat_mat, arma::cube& sxHat_cube_kh, arma::mat& Y0Hat_mat, arma::cube& Y0Hat_trunc_hg_cube, arma::mat& YTilde_mat, arma::mat& betaTilde_mat, arma::mat& UTilde_mat, arma::mat& Y0Tilde_mat, arma::mat& stat_nonstd_mat, arma::cube& stat_Sstd_cube_kh, arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec, int M_bts, double ap);
RcppExport SEXP _BTSinFLRM_PBouter(SEXP XSEXP, SEXP YSEXP, SEXP X0SEXP, SEXP tGridSEXP, SEXP betaHat_matSEXP, SEXP UHat_matSEXP, SEXP sxHat_cube_khSEXP, SEXP Y0Hat_matSEXP, SEXP Y0Hat_trunc_hg_cubeSEXP, SEXP YTilde_matSEXP, SEXP betaTilde_matSEXP, SEXP UTilde_matSEXP, SEXP Y0Tilde_matSEXP, SEXP stat_nonstd_matSEXP, SEXP stat_Sstd_cube_khSEXP, SEXP kn_vecSEXP, SEXP hn_vecSEXP, SEXP gn_vecSEXP, SEXP M_btsSEXP, SEXP apSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tGrid(tGridSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type betaHat_mat(betaHat_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type UHat_mat(UHat_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type sxHat_cube_kh(sxHat_cube_khSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y0Hat_mat(Y0Hat_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Y0Hat_trunc_hg_cube(Y0Hat_trunc_hg_cubeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type YTilde_mat(YTilde_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type betaTilde_mat(betaTilde_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type UTilde_mat(UTilde_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y0Tilde_mat(Y0Tilde_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type stat_nonstd_mat(stat_nonstd_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type stat_Sstd_cube_kh(stat_Sstd_cube_khSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type kn_vec(kn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type hn_vec(hn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gn_vec(gn_vecSEXP);
    Rcpp::traits::input_parameter< int >::type M_bts(M_btsSEXP);
    Rcpp::traits::input_parameter< double >::type ap(apSEXP);
    rcpp_result_gen = Rcpp::wrap(PBouter(X, Y, X0, tGrid, betaHat_mat, UHat_mat, sxHat_cube_kh, Y0Hat_mat, Y0Hat_trunc_hg_cube, YTilde_mat, betaTilde_mat, UTilde_mat, Y0Tilde_mat, stat_nonstd_mat, stat_Sstd_cube_kh, kn_vec, hn_vec, gn_vec, M_bts, ap));
    return rcpp_result_gen;
END_RCPP
}
// WBinner
List WBinner(arma::vec& W, arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, arma::mat& betaHat_mat, arma::mat& Y0Hat_mat, arma::cube& Y0Hat_trunc_hg_cube, arma::mat& erHat_mat, arma::mat& YHat_mat, arma::mat& XCent, arma::mat& YTilde_mat, arma::mat& Y0Tilde_mat, arma::cube& GaInvHat_cube, arma::cube& sxHat_cube_kh, arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec);
RcppExport SEXP _BTSinFLRM_WBinner(SEXP WSEXP, SEXP XSEXP, SEXP YSEXP, SEXP X0SEXP, SEXP tGridSEXP, SEXP betaHat_matSEXP, SEXP Y0Hat_matSEXP, SEXP Y0Hat_trunc_hg_cubeSEXP, SEXP erHat_matSEXP, SEXP YHat_matSEXP, SEXP XCentSEXP, SEXP YTilde_matSEXP, SEXP Y0Tilde_matSEXP, SEXP GaInvHat_cubeSEXP, SEXP sxHat_cube_khSEXP, SEXP kn_vecSEXP, SEXP hn_vecSEXP, SEXP gn_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tGrid(tGridSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type betaHat_mat(betaHat_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y0Hat_mat(Y0Hat_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Y0Hat_trunc_hg_cube(Y0Hat_trunc_hg_cubeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type erHat_mat(erHat_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type YHat_mat(YHat_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type XCent(XCentSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type YTilde_mat(YTilde_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y0Tilde_mat(Y0Tilde_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type GaInvHat_cube(GaInvHat_cubeSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type sxHat_cube_kh(sxHat_cube_khSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type kn_vec(kn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type hn_vec(hn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gn_vec(gn_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(WBinner(W, X, Y, X0, tGrid, betaHat_mat, Y0Hat_mat, Y0Hat_trunc_hg_cube, erHat_mat, YHat_mat, XCent, YTilde_mat, Y0Tilde_mat, GaInvHat_cube, sxHat_cube_kh, kn_vec, hn_vec, gn_vec));
    return rcpp_result_gen;
END_RCPP
}
// WBouter
List WBouter(arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, arma::mat& betaHat_mat, arma::mat& Y0Hat_mat, arma::cube& Y0Hat_trunc_hg_cube, arma::mat& erHat_mat, arma::mat& YHat_mat, arma::mat& XCent, arma::mat& YTilde_mat, arma::mat& Y0Tilde_mat, arma::cube& GaInvHat_cube, arma::cube& sxHat_cube_kh, arma::mat& stat_nonstd_mat, arma::cube& stat_Sstd_cube_kh, arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec, arma::mat& W, double ap);
RcppExport SEXP _BTSinFLRM_WBouter(SEXP XSEXP, SEXP YSEXP, SEXP X0SEXP, SEXP tGridSEXP, SEXP betaHat_matSEXP, SEXP Y0Hat_matSEXP, SEXP Y0Hat_trunc_hg_cubeSEXP, SEXP erHat_matSEXP, SEXP YHat_matSEXP, SEXP XCentSEXP, SEXP YTilde_matSEXP, SEXP Y0Tilde_matSEXP, SEXP GaInvHat_cubeSEXP, SEXP sxHat_cube_khSEXP, SEXP stat_nonstd_matSEXP, SEXP stat_Sstd_cube_khSEXP, SEXP kn_vecSEXP, SEXP hn_vecSEXP, SEXP gn_vecSEXP, SEXP WSEXP, SEXP apSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tGrid(tGridSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type betaHat_mat(betaHat_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y0Hat_mat(Y0Hat_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Y0Hat_trunc_hg_cube(Y0Hat_trunc_hg_cubeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type erHat_mat(erHat_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type YHat_mat(YHat_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type XCent(XCentSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type YTilde_mat(YTilde_matSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y0Tilde_mat(Y0Tilde_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type GaInvHat_cube(GaInvHat_cubeSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type sxHat_cube_kh(sxHat_cube_khSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type stat_nonstd_mat(stat_nonstd_matSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type stat_Sstd_cube_kh(stat_Sstd_cube_khSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type kn_vec(kn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type hn_vec(hn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gn_vec(gn_vecSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type ap(apSEXP);
    rcpp_result_gen = Rcpp::wrap(WBouter(X, Y, X0, tGrid, betaHat_mat, Y0Hat_mat, Y0Hat_trunc_hg_cube, erHat_mat, YHat_mat, XCent, YTilde_mat, Y0Tilde_mat, GaInvHat_cube, sxHat_cube_kh, stat_nonstd_mat, stat_Sstd_cube_kh, kn_vec, hn_vec, gn_vec, W, ap));
    return rcpp_result_gen;
END_RCPP
}