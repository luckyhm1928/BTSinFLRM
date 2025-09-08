// source_BST_FLR_paired.cpp


#include <RcppArmadillo.h>

#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

#define RCPP_ARMADILLO_FIX_Field

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat orthoL2Equidense(arma::mat& X, arma::vec& tGrid){
  int n = X.n_rows;
  int tt = X.n_cols;
  double scal = (range(tGrid) + tGrid(1) - tGrid(0)) / tt;  // left 
  
  arma::mat Xon(n, tt, fill::zeros);
  
  arma::vec x1 = trans(X.row(0));
  Xon.row(0) = trans(x1/sqrt(sum(square(x1))*scal));
  
  for(int i=1; i<n; ++i){
    arma::vec x_i = trans(X.row(i));
    arma::mat X_i1 = Xon.rows(0, i-1);
    arma::vec v_i = x_i - trans(X_i1) * (X_i1 * x_i * scal);
    Xon.row(i) = trans(v_i/sqrt(sum(square(v_i))*scal));
    
  }
  return(Xon);
  
}


// =========================================================
// simulation

//' @export
//[[Rcpp::export]]
arma::mat Xsim(int n, arma::vec& ev, arma::mat& ef, int xi_type=1, double nu_xi=1){
  int J = ev.size(); int tt = ef.n_cols;
  arma::mat X(n, tt, fill::zeros);
  if(xi_type==1){   // normal
    for(int i=0; i<n; i++){
      rowvec xi_vec = R::rnorm(0, 1) * Rcpp::rnorm(J, 0, 1); 
      rowvec X_i(tt, fill::zeros);
      for(int j=0; j<J; j++){
        X_i += sqrt(ev(j))*ef.row(j)*xi_vec(j);
      }
      X.row(i) = X_i;
    }
  }else if(xi_type==2){   // t distribution
    double sd = sqrt(nu_xi/(nu_xi-2));
    for(int i=0; i<n; i++){
      rowvec xi_vec =(R::rt(nu_xi)/sd) * Rcpp::rnorm(J, 0, 1); 
      rowvec X_i(tt, fill::zeros);
      for(int j=0; j<J; j++){
        X_i += sqrt(ev(j))*ef.row(j)*xi_vec(j);
      }
      X.row(i) = X_i;
    }
  }
  
  return(X);
}

//' @export
//[[Rcpp::export]]
List Ysim(arma::mat& X, arma::vec& beta, arma::vec& ev, arma::mat& ef, arma::vec& tGrid, 
          int distn_type=1, int var_type=1, int dep_type=1){
  int n = X.n_rows; int tt = X.n_cols;
  double scal = (range(tGrid) + tGrid(1) - tGrid(0)) / tt;
  
  arma::vec er(n, fill::zeros);
  if(distn_type==1){   // chi-square distribution
    if(var_type==1){   // cond'n var = norm
      if(dep_type==1){   // independent
        double nu = sum(ev)/2;
        er = Rcpp::rchisq(n, nu) - nu;
      }else if(dep_type==2){   // dependent
        arma::vec nu = sum(pow(X, 2), 1) * scal / 2;
        for(int i=0; i<n; ++i){
          er(i) = R::rchisq(nu(i)) - nu(i);
        }
      }
    }else if(var_type==2){   // cond'n var = inner product
      arma::vec rho(tt, fill::zeros);
      for(int t=0; t<tt; ++t){ rho(t) = pow(tGrid(t), 3) - 1.5 * tGrid(t) - 2.5; }
      
      if(dep_type==1){   // independent
        double nu = sum(ev % pow(ef * rho,2)) * scal * scal / 2;
        er = Rcpp::rchisq(n, nu) - nu;
      }else if(dep_type==2){   // dependent
        arma::vec nu = pow((X * rho), 2) * scal * scal / 2;
        for(int i=0; i<n; ++i){
          er(i) = R::rchisq(nu(i)) - nu(i);
        }
      }
    }
  }else if(distn_type==2){   // uniform distribution
    if(var_type==1){   // cond'n var = norm
      if(dep_type==1){   // independent
        double nu = sqrt(3*sum(ev));
        er = Rcpp::runif(n, -nu, nu);
      }else if(dep_type==2){   // dependent
        arma::vec nu = sqrt(3*sum(pow(X, 2), 1) * scal);
        for(int i=0; i<n; ++i){
          er(i) = R::runif(-nu(i), nu(i));
        }
      }
    }else if(var_type==2){   // cond'n var = inner product
      arma::vec rho(tt, fill::zeros);
      for(int t=0; t<tt; ++t){ rho(t) = pow(tGrid(t), 3) - 1.5 * tGrid(t) - 2.5; }
      
      if(dep_type==1){   // independent
        double nu = sqrt(3 * sum(ev % pow(ef * rho,2)) * scal * scal);
        er = Rcpp::runif(n, -nu, nu);
      }else if(dep_type==2){   // dependent
        arma::vec nu = sqrt(3 * pow((X * rho), 2) * scal * scal);
        for(int i=0; i<n; ++i){
          er(i) = R::runif(-nu(i), nu(i));
        }
      }
    }
  }
  
  arma::vec Xbeta = (X * beta) * scal;
  arma::vec Y = Xbeta + er;
  List res=List::create(Named("Xbeta") = Xbeta, Named("Y") = Y);
  return(res);
}






//[[Rcpp::export]]
List inferFLRM(arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, 
               arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec,
               double ap=0.05){
  
  arma::vec ap_vec = {1-ap};
  arma::vec qp_vec = {ap*0.5, 1-ap*0.5};
  double zval = R::qnorm(qp_vec(1), 0, 1, 1, 0);
  
  int n = X.n_rows; int n0 = X0.n_rows; int tt = X.n_cols;
  double scal = range(tGrid) / tt;
  int k_size = kn_vec.size(); 
  int h_size = hn_vec.size(); 
  int g_size = gn_vec.size(); 
  
  // arma::vec qp_vec = {ap*0.5, 1-ap*0.5};
  // arma::vec ap_vec = {1-ap};
  
  arma::vec tune(3); tune(0) = max(kn_vec); tune(1) = max(hn_vec); tune(2) = max(gn_vec);
  int tmx = tune.max();
  
  // orthonormalization of X0
  
  arma::mat X0on(n, tt, fill::zeros);
  
  arma::vec x01 = trans(X0.row(0));
  X0on.row(0) = trans(x01/sqrt(sum(square(x01))*scal));
  
  for(int i0=1; i0<n0; ++i0){
    arma::vec x_i0 = trans(X0.row(i0));
    arma::mat X_i01 = X0on.rows(0, i0-1);
    arma::vec v_i0 = x_i0 - trans(X_i01) * (X_i01 * x_i0 * scal);
    X0on.row(i0) = trans(v_i0/sqrt(sum(square(v_i0))*scal));
  }
  arma::mat projX0 = trans(X0on) * X0on;
  
  // centering
  arma::vec Xbar = trans(mean(X)); double Ybar = mean(Y);
  arma::vec YCent = Y - Ybar;
  arma::mat XCent(n,tt);
  arma::vec DtHat(tt);
  for(int t=0;t<tt;++t){
    arma::vec XCent_inner = X.col(t) - Xbar(t);
    XCent.col(t) = XCent_inner;
    DtHat(t) = mean(XCent_inner % YCent);
  }
  
  // estimation
  
  arma::mat GaHat = cov(XCent); arma::vec gaHat; arma::mat phiHat; 
  eig_sym(gaHat, phiHat, GaHat);
  gaHat = reverse(gaHat) * scal; 
  phiHat = fliplr(phiHat) / sqrt(scal);
  
  arma::cube PiHat_cube(tt, tt, tmx, fill::zeros);
  arma::cube GaInvHat_cube(tt, tt, tmx, fill::zeros);
  arma::mat temp(tt, tt, fill::zeros);
  arma::mat tempp(tt, tt, fill::zeros);
  for(int j=0; j<tmx; ++j){
    arma::mat piHat_j = phiHat.col(j) * trans(phiHat.col(j));
    temp = temp + piHat_j;
    PiHat_cube.slice(j) = temp;
    tempp = tempp + piHat_j / gaHat(j);
    GaInvHat_cube.slice(j) = tempp;
  }
  
  arma::mat betaHat_mat(tmx, tt, fill::zeros);
  arma::mat Y0Hat_mat(tmx, n0, fill::zeros);
  arma::mat erHat_mat(tmx, n, fill::zeros);
  arma::mat erHatCent_mat(tmx, n, fill::zeros);
  arma::mat YHat_mat(tmx, n, fill::zeros);
  arma::mat UHat_mat(tmx, tt, fill::zeros);
  arma::cube LdHat_cube(tt, tt, tmx, fill::zeros);
  
  arma::mat betaTilde_mat(tmx, tt, fill::zeros);
  arma::mat Y0Tilde_mat(tmx, n0, fill::zeros);
  arma::mat YTilde_mat(tmx, n, fill::zeros);
  arma::mat erTilde_mat(tmx, n, fill::zeros);
  arma::mat erTildeCent_mat(tmx, n, fill::zeros);
  arma::mat UTilde_mat(tmx, tt, fill::zeros);
  
  for(int j=0; j<tmx; ++j){
    arma::vec betaHat = GaInvHat_cube.slice(j) * DtHat * scal;
    betaHat_mat.row(j) = trans(betaHat);
    arma::vec Y0Hat_j = X0 * betaHat * scal;
    Y0Hat_mat.row(j) = trans(Y0Hat_j);
    
    arma::vec erHat = Y - X * betaHat * scal;
    arma::vec erHatCent = erHat - mean(erHat);
    erHat_mat.row(j) = trans(erHat);
    erHatCent_mat.row(j) = trans(erHatCent);
    
    // for homoscedasticity
    
    arma::vec YHat_j = X * betaHat * scal;
    YHat_mat.row(j) = trans(YHat_j);
    
    // for heteroscedasticity
    
    arma::vec UHat(tt, fill::zeros);
    for(int t=0;t<tt;++t){
      UHat(t) = mean(XCent.col(t) % erHatCent);
    }
    UHat_mat.row(j) = trans(UHat);
    
    arma::mat XCent_erHatCent(n, tt, fill::zeros);
    for(int i=0; i<n; ++i){
      XCent_erHatCent.row(i) = XCent.row(i) * erHatCent(i);
    }
    LdHat_cube.slice(j) = cov(XCent_erHatCent);
    
    // enforcing the null
    
    arma::vec PiX0betaHat = projX0 * betaHat * scal;
    // arma::vec PiX0betaHat = trans(X0) * Y0Hat_j;
    arma::vec betaTilde = betaHat-PiX0betaHat;
    betaTilde_mat.row(j) = trans(betaTilde);
    Y0Tilde_mat.row(j) = trans(X0 * betaTilde * scal);
    arma::vec YTilde = Y - X * PiX0betaHat * scal;
    YTilde_mat.row(j) = trans(YTilde);
    
    arma::vec erTilde = YTilde - X * betaTilde * scal;
    arma::vec erTildeCent = erTilde - mean(erTilde);
    erTilde_mat.row(j) = trans(erTilde);
    erTildeCent_mat.row(j) = trans(erTildeCent);
    
    arma::vec UTilde(tt, fill::zeros);
    for(int t=0;t<tt;++t){
      UTilde(t) = mean(XCent.col(t) % erTildeCent);
    }
    UTilde_mat.row(j) = trans(UTilde);
    
  }
  
  // for biased centering in the bootstrap world
  arma::cube Y0Hat_trunc_hg_cube(h_size, g_size, n0, fill::zeros);
  for(int ii_g=0; ii_g<g_size; ++ii_g){
    int g = gn_vec(ii_g);
    
    arma::vec betaHat_g = trans(betaHat_mat.row(g-1));
    
    for(int ii_h=0; ii_h<h_size; ++ii_h){
      int h = hn_vec(ii_h);
      
      arma::mat PiHat_h = PiHat_cube.slice(h-1);
      
      arma::vec Y0Hat_g_trunc_h = X0 * (PiHat_h * betaHat_g * scal) * scal;
      Y0Hat_trunc_hg_cube.tube(ii_h, ii_g) = Y0Hat_g_trunc_h;
    }
  }
  
  // results 
  
  List res_est = List::create(
    Named("GaInvHat") = GaInvHat_cube,
    Named("betaHat") = betaHat_mat,
    Named("Y0Hat") = Y0Hat_mat,
    Named("YHat") = YHat_mat,
    Named("erHat") = erHat_mat,
    Named("erHatCent") = erHatCent_mat,
    Named("UHat") = UHat_mat,
    Named("LdHat") = LdHat_cube,
    Named("Y0Hat_trunc_hg") = Y0Hat_trunc_hg_cube
  );
  
  List res_estTilde = List::create(
    Named("betaTilde") = betaTilde_mat,
    Named("Y0Tilde") = Y0Tilde_mat,
    Named("YTilde") = YTilde_mat,
    Named("erTilde") = erTilde_mat,
    Named("erTildeCent") = erTildeCent_mat,
    Named("UTilde") = UTilde_mat
  );
  
  
  // inference
  
  arma::mat stat_nonstd_mat_h(h_size, 2, fill::zeros);
  
  // homoscedasticity
  
  arma::vec sHatsq_vec_k(k_size, fill::zeros);
  arma::mat txHat_mat_h(h_size, n0, fill::zeros);
  arma::mat Tstd_mat_h(h_size, n0, fill::zeros);
  arma::mat stat_Tstd_mat_h(h_size, 2, fill::zeros);
  
  arma::cube upperCI_CLT_homo_cube_hk(k_size, h_size, n0, fill::zeros);
  arma::cube lowerCI_CLT_homo_cube_hk(k_size, h_size, n0, fill::zeros);
  arma::cube upperPI_CLT_homo_cube_hk(k_size, h_size, n0, fill::zeros);
  arma::cube lowerPI_CLT_homo_cube_hk(k_size, h_size, n0, fill::zeros);
  
  // heteroscedasticity
  
  arma::cube sxHat_cube_kh(k_size, h_size, n0, fill::zeros);
  arma::cube Sstd_cube_kh(k_size, h_size, n0, fill::zeros);
  arma::cube stat_Sstd_cube_kh(k_size, h_size, 2, fill::zeros);
  
  arma::cube upperCI_CLT_hetero_cube_hk(k_size, h_size, n0, fill::zeros);
  arma::cube lowerCI_CLT_hetero_cube_hk(k_size, h_size, n0, fill::zeros);
  
  for(int ii_h=0; ii_h<h_size; ++ii_h){
    int h = hn_vec(ii_h);
    
    // without studentization
    arma::vec Y0Hat_h = trans(Y0Hat_mat.row(h-1));
    stat_nonstd_mat_h(ii_h, 0) = sum(square(Y0Hat_h));
    stat_nonstd_mat_h(ii_h, 1) = max(abs(Y0Hat_h));
    
    // homoscedasticity
    arma::vec txHat_h(n0, fill::zeros);
    for(int j=0; j<h; j++){
      txHat_h += (X0 * phiHat.col(j)) % 
        (X0 * phiHat.col(j)) * pow(scal, 2) / gaHat(j);
    }
    txHat_mat_h.row(ii_h) = trans(txHat_h);
    
    // for testing
    arma::vec Tstd = sqrt(n) * Y0Hat_h / sqrt(txHat_h);
    Tstd_mat_h.row(ii_h) = trans(Tstd);
    stat_Tstd_mat_h(ii_h, 0) = sum(square(Tstd));
    stat_Tstd_mat_h(ii_h, 1) = max(abs(Tstd));
    
    
    // heteroscedasticity
    
    arma::mat GaInvHat_h = GaInvHat_cube.slice(h-1);
    arma::mat gix = GaInvHat_h * trans(X0) * scal;
    
    for(int ii_k=0; ii_k<k_size; ++ii_k){
      int k = kn_vec(ii_k);
      
      // homoscedasticity
      arma::vec erHatCent = trans(erHatCent_mat.row(k-1));
      double sHatsq = as_scalar(mean(pow(erHatCent,2)));
      
      // intervals
      
      arma::vec a_CI_homo = sqrt( sHatsq * txHat_h / n);
      arma::vec a_PI_homo = sqrt( sHatsq * (txHat_h / n + 1) );
      
      upperCI_CLT_homo_cube_hk.tube(ii_k,ii_h) = Y0Hat_h + zval * a_CI_homo;
      lowerCI_CLT_homo_cube_hk.tube(ii_k,ii_h) = Y0Hat_h - zval * a_CI_homo;
      upperPI_CLT_homo_cube_hk.tube(ii_k,ii_h) = Y0Hat_h + zval * a_PI_homo;
      lowerPI_CLT_homo_cube_hk.tube(ii_k,ii_h) = Y0Hat_h + zval * a_PI_homo;
      
      // heteroscedasticity
      
      arma::vec sxHat_h = trans(sum((LdHat_cube.slice(k-1) * gix * scal) % gix * scal));
      sxHat_cube_kh.tube(ii_k, ii_h) = sxHat_h;
      
      arma::vec a_CI_hetero = sqrt(sxHat_h / n);
      upperCI_CLT_hetero_cube_hk.tube(ii_k,ii_h) = Y0Hat_h + zval * a_CI_hetero;
      lowerCI_CLT_hetero_cube_hk.tube(ii_k,ii_h) = Y0Hat_h - zval * a_CI_hetero;
      // no PI by CLT under hetero
      
      // for testing
      arma::vec Sstd = sqrt(n) * Y0Hat_h / sqrt(sxHat_h);
      
      Sstd_cube_kh.tube(ii_k, ii_h) = Sstd;
      stat_Sstd_cube_kh(ii_k, ii_h, 0) = sum(square(Sstd));
      stat_Sstd_cube_kh(ii_k, ii_h, 1) = max(abs(Sstd));
    }
  }
  
  
  // results
  
  List resCLThomo = List::create(
    Named("upperCI_CLT_homo") = upperCI_CLT_homo_cube_hk,
    Named("lowerCI_CLT_homo") = lowerCI_CLT_homo_cube_hk,
    Named("upperPI_CLT_homo") = upperPI_CLT_homo_cube_hk,
    Named("lowerPI_CLT_homo") = lowerPI_CLT_homo_cube_hk,
    Named("sHatsq") = sHatsq_vec_k,
    Named("txHat") = txHat_mat_h
  );
  
  List resCLThetero = List::create(
    Named("upperCI_CLT_hetero") = upperCI_CLT_hetero_cube_hk,
    Named("lowerCI_CLT_hetero") = lowerCI_CLT_hetero_cube_hk,
    Named("sxHat") = sxHat_cube_kh
  );
  
  List resSTAThomo = List::create(
    Named("stat_nonstd") = stat_nonstd_mat_h,
    Named("Tstd") = Tstd_mat_h,
    Named("stat_Tstd") = stat_Tstd_mat_h
  );
  
  List resSTAThetero = List::create(
    Named("stat_nonstd") = stat_nonstd_mat_h,
    Named("Sstd") = Sstd_cube_kh,
    Named("stat_Sstd") = stat_Sstd_cube_kh
  );
  
  
  
  // final results
  
  List res = List::create(
    Named("est") = res_est,
    Named("estTilde") = res_estTilde,
    Named("CLThomo") = resCLThomo,
    Named("CLThetero") = resCLThetero,
    Named("STAThomo") = resSTAThomo,
    Named("STAThetero") = resSTAThetero
  );
  
  return(res);
}














//[[Rcpp::export]]
List RBinner(
    arma::vec& idx, arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, 
    arma::mat& betaHat_mat, arma::mat& Y0Hat_mat, arma::cube& Y0Hat_trunc_hg_cube,
    arma::mat& erHatCent_mat, arma::mat& YHat_mat, arma::mat& XCent, 
    arma::mat& YTilde_mat, arma::mat& Y0Tilde_mat, arma::cube& GaInvHat_cube, arma::mat& txHat_mat_h,
    arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec
){
  
  int n = X.n_rows; int n0 = X0.n_rows; int tt = X.n_cols;
  double scal = range(tGrid) / tt;
  int k_size = kn_vec.size(); int h_size = hn_vec.size(); int g_size = gn_vec.size(); 
  arma::vec tune(3); tune(0) = max(kn_vec); tune(1) = max(hn_vec); tune(2) = max(gn_vec);
  int tmx = tune.max();
  
  // residual bootstrap
  
  arma::mat erStar_mat(k_size, n, fill::zeros);
  arma::mat erStarCent_mat(k_size, n, fill::zeros);
  
  arma::cube YStar_cube(k_size, g_size, n, fill::zeros);
  arma::cube YStarCent_cube(k_size, g_size, n, fill::zeros);
  
  arma::cube Y0Star_cube(k_size, g_size, n0, fill::zeros);
  
  // arma::field<arma::vec> Y0HatStar_field(k_size, h_size, g_size);
  
  arma::field<arma::vec> rootsCI_trunc_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsCI_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsPI_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsStdCI_trunc_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsStdCI_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsStdPI_field_inner(k_size, h_size, g_size);
  
  arma::field<arma::vec> statStar_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> statStar_std_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> statTildeStar_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> statTildeStar_std_field_inner(k_size, h_size, g_size);
  
  for(int ii_k=0; ii_k<k_size; ++ii_k){
    int k = kn_vec(ii_k);
    
    arma::vec erHatCent = trans(erHatCent_mat.row(k-1));
    arma::vec erStar(n);
    for(int i=0; i<n; ++i){ erStar(i) = erHatCent(idx(i)); }
    arma::vec erStarCent = erStar - mean(erStar);
    
    erStar_mat.row(ii_k) = trans(erStar);
    erStarCent_mat.row(ii_k) = trans(erStarCent);
    
    // prediction
    arma::vec idx0 = floor(Rcpp::runif(n0)*n);
    arma::vec er0Star(n0);
    for(int i0=0; i0<n0; ++i0){ er0Star(i0) = erHatCent(idx0(i0)); }
    
    for(int ii_g=0; ii_g<g_size; ++ii_g){
      int g = gn_vec(ii_g);
      
      arma::vec YHat_g = trans(YHat_mat.row(g-1));
      arma::vec YStar = YHat_g + erStar;
      double YStarbar = mean(YStar);
      arma::vec YStarCent = YStar - YStarbar;
      
      YStar_cube.tube(ii_k, ii_g) = YStar;
      YStarCent_cube.tube(ii_k, ii_g) = YStarCent;
      
      arma::vec DtHatStar(tt);
      for(int t=0;t<tt;++t){
        DtHatStar(t) = mean(XCent.col(t) % YStarCent);
      }
      
      // enforcing the null
      
      arma::vec YTilde_g = trans(YTilde_mat.row(g-1));
      arma::vec YTildeStar = YTilde_g + erStar;
      double YTildeStarbar = mean(YTildeStar);
      arma::vec YTildeStarCent = YTildeStar - YTildeStarbar;
      
      YStar_cube.tube(ii_k, ii_g) = YStar;
      YStarCent_cube.tube(ii_k, ii_g) = YStarCent;
      
      arma::vec DtTildeStar(tt);
      for(int t=0;t<tt;++t){
        DtTildeStar(t) = mean(XCent.col(t) % YTildeStarCent);
      }
      
      arma::vec Y0Tilde_g = trans(Y0Tilde_mat.row(g-1));
      
      // prediction
      arma::vec Y0Hat_g = trans(Y0Hat_mat.row(g-1));
      arma::vec Y0Star = Y0Hat_g + er0Star;
      Y0Star_cube.tube(ii_k, ii_g) = Y0Star;
      
      for(int ii_h=0; ii_h<h_size; ++ii_h){
        int h = hn_vec(ii_h);
        
        arma::mat GaInvHat_h = GaInvHat_cube.slice(h-1);
        arma::vec betaHatStar_h = GaInvHat_h * DtHatStar * scal;
        arma::vec Y0HatStar_h = X0 * betaHatStar_h * scal;
        
        // enforcing the null
        arma::vec betaTildeStar_h = GaInvHat_h * DtTildeStar * scal;
        arma::vec Y0TildeStar_h = X0 * betaTildeStar_h * scal;
        
        
        // Y0HatStar_field(ii_k, ii_h, ii_g) = Y0HatStar_field;
        
        // store for bootstrap intervals
        
        arma::vec Y0Hat_g_trunc_h = Y0Hat_trunc_hg_cube.tube(ii_h, ii_g);
        arma::vec rootsCI_trunc = Y0HatStar_h - Y0Hat_g_trunc_h;
        arma::vec TStar = Y0HatStar_h - Y0Hat_g;
        arma::vec rootsPI = Y0HatStar_h - Y0Star;
        
        rootsCI_trunc_field_inner(ii_k, ii_h, ii_g) = rootsCI_trunc;
        rootsCI_field_inner(ii_k, ii_h, ii_g) = TStar;
        rootsPI_field_inner(ii_k, ii_h, ii_g) = rootsPI;
        
        arma::vec txHat_h = trans(txHat_mat_h.row(ii_h));
        rootsStdCI_trunc_field_inner(ii_k, ii_h, ii_g) = 
          sqrt(n) * rootsCI_trunc / sqrt(txHat_h);
        arma::vec TStar_std = sqrt(n) * TStar / sqrt(txHat_h);
        rootsStdCI_field_inner(ii_k, ii_h, ii_g) = TStar_std;
        rootsStdPI_field_inner(ii_k, ii_h, ii_g) = 
          sqrt(n) * rootsPI / sqrt(txHat_h/n+1);
        
        // testing
        
        // not enforcing the null
        
        double statL2Star = sum(square(TStar));
        double statMaxStar = max(abs(TStar));
        statStar_field_inner(ii_k, ii_h, ii_g) = {statL2Star, statMaxStar};
        
        double statL2Star_std = sum(square(TStar_std));
        double statMaxStar_std = max(abs(TStar_std));
        statStar_std_field_inner(ii_k, ii_h, ii_g) = {statL2Star_std, statMaxStar_std};
        
        // enforcing the null
        arma::vec TTildeStar = Y0TildeStar_h - Y0Tilde_g;
        double statTildeL2Star = sum(square(TTildeStar));
        double statTildeMaxStar = max(abs(TTildeStar));
        statTildeStar_field_inner(ii_k, ii_h, ii_g) = 
          {statTildeL2Star, statTildeMaxStar};
        
        arma::vec TTildeStar_std = sqrt(n) * TTildeStar / sqrt(txHat_h);
        double statTildeL2Star_std = sum(square(TTildeStar_std));
        double statTildeMaxStar_std = max(abs(TTildeStar_std));
        statTildeStar_std_field_inner(ii_k, ii_h, ii_g) = 
          {statTildeL2Star_std, statTildeMaxStar_std};
        
        
        
      }
    }
  }
  
  List res_rb_inner = List::create(
    Named("rootsCI_trunc_field_inner") = rootsCI_trunc_field_inner,
    Named("rootsCI_field_inner") = rootsCI_field_inner,
    Named("rootsPI_field_inner") = rootsPI_field_inner,
    Named("rootsStdCI_trunc_field_inner") = rootsStdCI_trunc_field_inner,
    Named("rootsStdCI_field_inner") = rootsStdCI_field_inner,
    Named("rootsStdPI_field_inner") = rootsStdPI_field_inner,
    Named("statStar_field_inner") = statStar_field_inner,
    Named("statStar_std_field_inner") = statStar_std_field_inner,
    Named("statTildeStar_field_inner") = statTildeStar_field_inner,
    Named("statTildeStar_std_field_inner") = statTildeStar_std_field_inner
  );
  
  return(res_rb_inner);
  
}





//[[Rcpp::export]]
List RBouter(
    arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, 
    arma::mat& betaHat_mat, arma::mat& Y0Hat_mat, arma::cube& Y0Hat_trunc_hg_cube,
    arma::mat& erHatCent_mat, arma::mat& YHat_mat, arma::mat& XCent, 
    arma::mat& YTilde_mat, arma::mat& Y0Tilde_mat, arma::cube& GaInvHat_cube, arma::mat& txHat_mat_h,
    arma::mat& stat_nonstd_mat, arma::mat& stat_Tstd_mat,
    arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec,
    int M_bts, double ap = 0.05){
  
  arma::vec ap_vec = {1-ap};
  arma::vec qp_vec = {ap*0.5, 1-ap*0.5};
  
  int n = X.n_rows; int n0 = X0.n_rows; int tt = X.n_cols;
  double scal = range(tGrid) / tt;
  int k_size = kn_vec.size(); 
  int h_size = hn_vec.size(); 
  int g_size = gn_vec.size(); 
  
  // arma::vec qp_vec = {ap*0.5, 1-ap*0.5};
  // arma::vec ap_vec = {1-ap};
  
  arma::vec tune(3); tune(0) = max(kn_vec); tune(1) = max(hn_vec); tune(2) = max(gn_vec);
  int tmx = tune.max();
  
  // 
  
  arma::field<arma::field<arma::vec>> rootsCI_trunc_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsCI_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsPI_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsStdCI_trunc_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsStdCI_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsStdPI_field(M_bts);
  
  arma::field<arma::field<arma::vec>> statStar_field(M_bts);
  arma::field<arma::field<arma::vec>> statStar_std_field(M_bts);
  arma::field<arma::field<arma::vec>> statTildeStar_field(M_bts);
  arma::field<arma::field<arma::vec>> statTildeStar_std_field(M_bts);
  
  for(int m=0; m<M_bts; ++m){
    
    // resampling
    arma::vec idx = floor(Rcpp::runif(n)*n);
    
    List resRBinner = 
      RBinner(
        idx, X, Y, X0, tGrid, 
        betaHat_mat, Y0Hat_mat, Y0Hat_trunc_hg_cube,
        erHatCent_mat, YHat_mat, XCent, 
        YTilde_mat, Y0Tilde_mat, GaInvHat_cube, txHat_mat_h,
        kn_vec, hn_vec, gn_vec
      );
    
    
    
    arma::field<arma::vec> rootsCI_trunc_field_inner =
      resRBinner["rootsCI_trunc_field_inner"];
    arma::field<arma::vec> rootsCI_field_inner = 
      resRBinner["rootsCI_field_inner"];
    arma::field<arma::vec> rootsPI_field_inner = 
      resRBinner["rootsPI_field_inner"];
    arma::field<arma::vec> rootsStdCI_trunc_field_inner = 
      resRBinner["rootsStdCI_trunc_field_inner"];
    arma::field<arma::vec> rootsStdCI_field_inner = 
      resRBinner["rootsStdCI_field_inner"];
    arma::field<arma::vec> rootsStdPI_field_inner = 
      resRBinner["rootsStdPI_field_inner"];
    
    
    rootsCI_trunc_field(m) = rootsCI_trunc_field_inner;
    rootsCI_field(m) = rootsCI_field_inner;
    rootsPI_field(m) = rootsPI_field_inner;
    rootsStdCI_trunc_field(m) = rootsStdCI_trunc_field_inner;
    rootsStdCI_field(m) = rootsStdCI_field_inner;
    rootsStdPI_field(m) = rootsStdPI_field_inner;
    
    arma::field<arma::vec> statStar_field_inner = 
      resRBinner["statStar_field_inner"];
    arma::field<arma::vec> statStar_std_field_inner = 
      resRBinner["statStar_std_field_inner"];
    arma::field<arma::vec> statTildeStar_field_inner = 
      resRBinner["statTildeStar_field_inner"];
    arma::field<arma::vec> statTildeStar_std_field_inner = 
      resRBinner["statTildeStar_std_field_inner"];
    
    statStar_field(m) = statStar_field_inner;
    statStar_std_field(m) = statStar_std_field_inner;
    statTildeStar_field(m) = statTildeStar_field_inner;
    statTildeStar_std_field(m) = statTildeStar_std_field_inner;
    
  }
  
  arma::field<arma::cube> roots_field_khg(k_size, h_size, g_size);
  arma::field<arma::cube> intervals_field_khg(k_size, h_size, g_size);
  arma::field<arma::cube> stat_field_khg(k_size, h_size, g_size);
  arma::field<arma::mat> pVal_field_khg(k_size, h_size, g_size);
  
  for(int ii_k=0; ii_k<k_size; ++ii_k){
    int k = kn_vec(ii_k);
    for(int ii_h=0; ii_h<h_size; ++ii_h){
      int h = hn_vec(ii_h);
      arma::vec Y0Hat_h = trans(Y0Hat_mat.row(h-1));
      arma::vec txHat_h = trans(txHat_mat_h.row(ii_h));
      
      for(int ii_g=0; ii_g<g_size; ++ii_g){
        int g = gn_vec(ii_g);
        
        arma::cube roots_cube(6, M_bts, n0, fill::zeros);
        // Hat_nonstd, Hat_std, Tilde_nonstd, Tilde_std
        arma::cube stat_cube(4, M_bts, 2, fill::zeros);
        arma::cube pVal_cube(4, M_bts, 2, fill::zeros);
        
        
        for(int m=0; m<M_bts; ++m){
          
          // intervals without std
          
          arma::field<arma::vec> rootsCI_trunc_field_inner = rootsCI_trunc_field(m);
          arma::vec rootsCI_trunc_vec = rootsCI_trunc_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(0, m) = rootsCI_trunc_vec;
          
          arma::field<arma::vec> rootsCI_field_inner = rootsCI_field(m);
          arma::vec rootsCI_vec = rootsCI_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(1, m) = rootsCI_vec;
          
          arma::field<arma::vec> rootsPI_field_inner = rootsPI_field(m);
          arma::vec rootsPI_vec = rootsPI_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(2, m) = rootsPI_vec;
          
          // intervals with std
          
          arma::field<arma::vec> rootsStdCI_trunc_field_inner = rootsStdCI_trunc_field(m);
          arma::vec rootsStdCI_trunc_vec = rootsStdCI_trunc_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(3, m) = rootsStdCI_trunc_vec;
          
          arma::field<arma::vec> rootsStdCI_field_inner = rootsStdCI_field(m);
          arma::vec rootsStdCI_vec = rootsStdCI_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(4, m) = rootsStdCI_vec;
          
          arma::field<arma::vec> rootsStdPI_field_inner = rootsStdPI_field(m);
          arma::vec rootsStdPI_vec = rootsStdPI_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(5, m) = rootsStdPI_vec;
          
          // testing: not enforcing the null
          
          arma::field<arma::vec> statStar_field_inner = statStar_field(m);
          arma::vec statStar_vec = statStar_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(0, m) = statStar_vec;
          
          arma::field<arma::vec> statStar_std_field_inner = statStar_std_field(m);
          arma::vec statStar_std_vec = statStar_std_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(1, m) = statStar_std_vec;
          
          if(stat_nonstd_mat(ii_h, 0)<statStar_vec(0)){pVal_cube(0, m, 0)=1;}
          if(stat_nonstd_mat(ii_h, 1)<statStar_vec(1)){pVal_cube(0, m, 1)=1;}
          if(stat_Tstd_mat(ii_h, 0)<statStar_std_vec(0)){pVal_cube(1, m, 0)=1;}
          if(stat_Tstd_mat(ii_h, 1)<statStar_std_vec(1)){pVal_cube(1, m, 1)=1;}
          
          // testing: enforcing the null
          
          arma::field<arma::vec> statTildeStar_field_inner = statTildeStar_field(m);
          arma::vec statTildeStar_vec = statTildeStar_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(2, m) = statTildeStar_vec;
          
          arma::field<arma::vec> statTildeStar_std_field_inner = statTildeStar_std_field(m);
          arma::vec statTildeStar_std_vec = statTildeStar_std_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(3, m) = statTildeStar_std_vec;
          
          if(stat_nonstd_mat(ii_h, 0)<statTildeStar_vec(0)){pVal_cube(2, m, 0)=1;}
          if(stat_nonstd_mat(ii_h, 1)<statTildeStar_vec(1)){pVal_cube(2, m, 0)=1;}
          if(stat_Tstd_mat(ii_h, 0)<statTildeStar_std_vec(0)){pVal_cube(3, m, 0)=1;}
          if(stat_Tstd_mat(ii_h, 1)<statTildeStar_std_vec(1)){pVal_cube(3, m, 1)=1;}
        }
        
        // ======
        // intervals
        
        roots_field_khg(ii_k, ii_h, ii_g) = roots_cube;
        arma::cube roots_cube_abs = abs(roots_cube);
        
        // individual intervals
        
        arma::cube II(n0, 2, 6);
        for(int l=0; l<3; ++l){
          // arma::mat temp = roots_cube.row(l);
          // arma::mat qrt = quantile(temp, qp_vec, 0);
          // arma::mat tempp = roots_cube_abs.row(l);
          // arma::mat qrt_sym = quantile(tempp, ap_vec, 0);
          for(int i0=0; i0<n0; ++i0){
            
            
            arma::mat temp = roots_cube.slice(i0);
            arma::rowvec tempp = temp.row(l);
            arma::rowvec qrt = quantile(tempp, qp_vec);
            
            arma::mat temp_abs = roots_cube_abs.slice(i0);
            arma::rowvec tempp_abs = temp_abs.row(l);
            arma::rowvec qrt_sym = quantile(tempp_abs, ap_vec);
            
            
            II(i0, 0, l) = Y0Hat_h(i0) - qrt(1);
            II(i0, 1, l) = Y0Hat_h(i0) - qrt(0);
            
            II(i0, 0, l+3) = Y0Hat_h(i0) - qrt_sym(0);
            II(i0, 1, l+3) = Y0Hat_h(i0) + qrt_sym(0);
          }
        }
        
        // simultaneous intervals
        
        arma::vec qrt_max_rb(6, fill::zeros);
        for(int l=0; l<6; ++l){
          arma::vec temp = max(roots_cube_abs.row(l), 1);
          qrt_max_rb(l) = as_scalar(quantile(temp, ap_vec));
        }
        
        arma::cube SI(n0, 2, 6);
        for(int i0=0; i0<n0; ++i0){
          for(int l=0; l<3; ++l){
            double qrt_max = qrt_max_rb(l);
            SI(i0, 0, l) = Y0Hat_h(i0) - qrt_max;
            SI(i0, 1, l) = Y0Hat_h(i0) + qrt_max;
          }
          double qrt_SI_std_rb_4 = qrt_max_rb(3) * sqrt(txHat_h(i0)/n);
          double qrt_SI_std_rb_5 = qrt_max_rb(4) * sqrt(txHat_h(i0)/n);
          double qrt_SI_std_rb_6 = qrt_max_rb(5) * sqrt(txHat_h(i0)/n+1);
          
          SI(i0, 0, 3) = Y0Hat_h(i0) - qrt_SI_std_rb_4;
          SI(i0, 1, 3) = Y0Hat_h(i0) + qrt_SI_std_rb_4;
          SI(i0, 0, 4) = Y0Hat_h(i0) - qrt_SI_std_rb_5;
          SI(i0, 1, 4) = Y0Hat_h(i0) + qrt_SI_std_rb_5;
          SI(i0, 0, 5) = Y0Hat_h(i0) - qrt_SI_std_rb_6;
          SI(i0, 1, 5) = Y0Hat_h(i0) + qrt_SI_std_rb_6;
          
        }
        
        // store
        
        intervals_field_khg(ii_k, ii_h, ii_g) = join_slices(II, SI);
        
        // ======
        // testing
        
        stat_field_khg(ii_k, ii_h, ii_g) = stat_cube;
        
        pVal_field_khg(ii_k, ii_h, ii_g) = mean(pVal_cube, 1);
        
        
      }
    }
  }
  
  List resRBouter = List::create(
    Named("roots") = roots_field_khg,
    Named("intervals") = intervals_field_khg,
    Named("stat") = stat_field_khg,
    Named("pVal") = pVal_field_khg
  );
  
  return(resRBouter);
  
}


















//[[Rcpp::export]]
List PBinner1(
    arma::vec& idx, arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, 
    arma::mat& betaHat_mat, arma::mat& UHat_mat, arma::cube& sxHat_cube_kh,
    arma::mat& Y0Hat_mat, arma::cube& Y0Hat_trunc_hg_cube,
    arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec
){
  
  int n = X.n_rows; int n0 = X0.n_rows; int tt = X.n_cols;
  double scal = range(tGrid) / tt;
  int k_size = kn_vec.size(); int h_size = hn_vec.size(); int g_size = gn_vec.size(); 
  arma::vec tune(3); tune(0) = max(kn_vec); tune(1) = max(hn_vec); tune(2) = max(gn_vec);
  int tmx = tune.max();
  
  arma::vec YStar(n);
  arma::mat XStar(n, tt);
  for(int i=0; i<n; ++i){
    YStar(i) = Y(idx(i));
    XStar.row(i) = X.row(idx(i));
  }
  
  
  arma::mat erStar_mat(tmx, n, fill::zeros);
  arma::mat erStarCent_mat(tmx, n, fill::zeros);
  
  for(int j=0; j<tmx; ++j){
    arma::vec betaHat = trans(betaHat_mat.row(j));
    
    arma::vec erStar = YStar - XStar * betaHat * scal;
    arma::vec erStarCent = erStar - mean(erStar);
    erStar_mat.row(j) = trans(erStar);
    erStarCent_mat.row(j) = trans(erStarCent);
  }
  
  // centering
  arma::vec XStarbar = trans(mean(XStar)); double YStarbar = mean(YStar);
  arma::vec YStarCent = YStar - YStarbar;
  arma::mat XStarCent(n,tt);
  arma::vec DtHatStar(tt);
  
  arma::mat X0Cent_by_star(n0,tt); // center X0 by \bar{X}^*
  for(int t=0;t<tt;++t){
    arma::vec XStarCent_inner = XStar.col(t) - XStarbar(t);
    XStarCent.col(t) = XStarCent_inner;
    DtHatStar(t) = mean(XStarCent_inner % YStarCent);
    
    arma::vec X0Cent_by_star_inner = X0.col(t) - XStarbar(t);
    X0Cent_by_star.col(t) = X0Cent_by_star_inner; // center X0 by \bar{X}^*
  }
  // X0 = X0Cent_by_star; // center X0 by \bar{X}^*
  
  // estimation
  
  arma::mat GaHatStar = cov(XStarCent); arma::vec gaHatStar; arma::mat phiHatStar; 
  eig_sym(gaHatStar, phiHatStar, GaHatStar);
  gaHatStar = reverse(gaHatStar) * scal; 
  phiHatStar = fliplr(phiHatStar) / sqrt(scal);
  
  arma::cube GaInvHatStar_cube(tt, tt, tmx, fill::zeros);
  arma::mat temppStar(tt, tt, fill::zeros);
  for(int j=0; j<tmx; ++j){
    temppStar = temppStar +
      phiHatStar.col(j) * trans(phiHatStar.col(j)) / gaHatStar(j);
    GaInvHatStar_cube.slice(j) = temppStar;
  }
  
  arma::cube betaHatStar_cube_g(tmx, g_size, tt, fill::zeros);
  arma::cube Y0HatStar_cube_g(tmx, g_size, n0, fill::zeros);
  arma::cube erHatStar_cube_g(tmx, g_size, n, fill::zeros);
  arma::cube erHatStarCent_cube_g(tmx, g_size, n, fill::zeros);
  arma::field<arma::mat> LdHatStar_field_g(tmx, g_size);
  
  
  for(int ii_g=0; ii_g<g_size; ++ii_g){
    int g = gn_vec(ii_g);
    arma::vec UHat = trans(UHat_mat.row(g-1));
    
    for(int j=0; j<tmx; ++j){
      arma::vec betaHatStar = GaInvHatStar_cube.slice(j) * (DtHatStar-UHat) * scal;
      betaHatStar_cube_g.tube(j, ii_g) = betaHatStar;
      arma::vec Y0HatStar_j = X0Cent_by_star * betaHatStar * scal; // X0 appear
      Y0HatStar_cube_g.tube(j, ii_g) = Y0HatStar_j;
      
      arma::vec erHatStar = YStar - XStar * betaHatStar * scal;
      arma::vec erHatStarCent = erHatStar - mean(erHatStar);
      erHatStar_cube_g.tube(j, ii_g) = trans(erHatStar);
      erHatStarCent_cube_g.tube(j, ii_g) = trans(erHatStarCent);
      
      arma::mat XStarCent_erHatStarCent(n, tt, fill::zeros);
      for(int i=0; i<n; ++i){
        XStarCent_erHatStarCent.row(i) = XStarCent.row(i) * erHatStarCent(i);
      }
      LdHatStar_field_g(j, ii_g) = cov(XStarCent_erHatStarCent);
      
    }
  }
  
  // inference
  
  // arma::vec idx0 = floor(Rcpp::runif(n0)*n);   // for prediction
  
  arma::field<arma::vec> rootsCI_trunc_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsCI_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsStd1CI_trunc_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsStd1CI_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsStd2CI_trunc_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsStd2CI_field_inner(k_size, h_size, g_size);
  
  arma::field<arma::vec> statStar_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> statStar_std1_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> statStar_std2_field_inner(k_size, h_size, g_size);
  
  for(int ii_k=0; ii_k<k_size; ++ii_k){
    int k = kn_vec(ii_k);
    for(int ii_g=0; ii_g<g_size; ++ii_g){
      int g = gn_vec(ii_g);
      arma::mat LdHatStar = LdHatStar_field_g(k-1, ii_g);
      arma::vec Y0Hat_g = trans(Y0Hat_mat.row(g-1));
      
      // arma::vec er0Star(n0);
      // for(int i0=0; i0<n0; ++i0){ er0Star(i0) = erStarCent(idx0(i0)); }
      // arma::vec Y0Star = Y0Hat_g + er0Star;
      
      for(int ii_h=0; ii_h<h_size; ++ii_h){
        int h = hn_vec(ii_h);
        
        arma::vec Y0HatStar_h = Y0HatStar_cube_g.tube(h-1, ii_g);
        
        
        arma::vec Y0Hat_g_trunc_h = Y0Hat_trunc_hg_cube.tube(ii_h, ii_g);
        arma::vec rootsCI_trunc = Y0HatStar_h - Y0Hat_g_trunc_h;
        arma::vec TStar = Y0HatStar_h - Y0Hat_g;
        // arma::vec rootsPI = Y0HatStar_h - Y0Star;
        
        rootsCI_trunc_field_inner(ii_k, ii_h, ii_g) = rootsCI_trunc;
        rootsCI_field_inner(ii_k, ii_h, ii_g) = TStar;
        // rootsPI_field_inner(ii_k, ii_h, ii_g) = rootsPI;
        
        
        
        // std 1
        arma::vec sxHat_h = sxHat_cube_kh.tube(ii_k, ii_h) ; // X0 appear, but centered by \bar{X}
        arma::vec TStar_std1 = sqrt(n) * TStar / sqrt(sxHat_h);
        rootsStd1CI_trunc_field_inner(ii_k, ii_h, ii_g) = 
          sqrt(n) * rootsCI_trunc / sqrt(sxHat_h);
        rootsStd1CI_field_inner(ii_k, ii_h, ii_g) = TStar_std1;
        
        
        // std 2
        
        arma::mat GaInvHatStar_h = GaInvHatStar_cube.slice(h-1);
        arma::mat gixStar = GaInvHatStar_h * trans(X0Cent_by_star) * scal; // X0 appear
        arma::vec sxHatStar_h = trans(sum((LdHatStar * gixStar * scal) % gixStar * scal));
        
        arma::vec TStar_std2 = sqrt(n) * TStar / sqrt(sxHatStar_h);
        rootsStd2CI_trunc_field_inner(ii_k, ii_h, ii_g) = 
          sqrt(n) * rootsCI_trunc / sqrt(sxHatStar_h);
        rootsStd2CI_field_inner(ii_k, ii_h, ii_g) = TStar_std2;
        
        
        
        
        // testing
        
        double statL2Star = sum(square(TStar));
        double statMaxStar = max(abs(TStar));
        statStar_field_inner(ii_k, ii_h, ii_g) =
          {statL2Star, statMaxStar};
        
        double statL2Star_std1 = sum(square(TStar_std1));
        double statMaxStar_std1 = max(abs(TStar_std1));
        statStar_std1_field_inner(ii_k, ii_h, ii_g) = 
          {statL2Star_std1, statMaxStar_std1};
        
        double statL2Star_std2 = sum(square(TStar_std2));
        double statMaxStar_std2 = max(abs(TStar_std2));
        statStar_std2_field_inner(ii_k, ii_h, ii_g) = 
          {statL2Star_std2, statMaxStar_std2};
        
        
      }
    }
  }
  
  List res_pb_inner = List::create(
    Named("rootsCI_trunc_field_inner") = rootsCI_trunc_field_inner,
    Named("rootsCI_field_inner") = rootsCI_field_inner,
    Named("rootsStd1CI_trunc_field_inner") = rootsStd1CI_trunc_field_inner,
    Named("rootsStd1CI_field_inner") = rootsStd1CI_field_inner,
    Named("rootsStd2CI_trunc_field_inner") = rootsStd2CI_trunc_field_inner,
    Named("rootsStd2CI_field_inner") = rootsStd2CI_field_inner,
    Named("statStar_field_inner") = statStar_field_inner,
    Named("statStar_std1_field_inner") = statStar_std1_field_inner,
    Named("statStar_std2_field_inner") = statStar_std2_field_inner
  );
  
  return(res_pb_inner);
}





//[[Rcpp::export]]
List PBinner2(
    arma::vec& idxTilde, arma::mat& X, arma::mat& YTilde_mat, arma::mat& X0, arma::vec& tGrid, 
    arma::mat& betaTilde_mat, arma::mat& UTilde_mat, arma::cube& sxHat_cube_kh,
    arma::mat& Y0Tilde_mat, 
    arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec
){
  
  int n = X.n_rows; int n0 = X0.n_rows; int tt = X.n_cols;
  double scal = range(tGrid) / tt;
  int k_size = kn_vec.size(); int h_size = hn_vec.size(); int g_size = gn_vec.size(); 
  arma::vec tune(3); tune(0) = max(kn_vec); tune(1) = max(hn_vec); tune(2) = max(gn_vec);
  int tmx = tune.max();
  
  arma::mat YTildeStar_mat(tmx, n, fill::zeros);
  arma::mat XTildeStar(n, tt);
  for(int i=0; i<n; ++i){
    YTildeStar_mat.col(i) = YTilde_mat.col(idxTilde(i));
    XTildeStar.row(i) = X.row(idxTilde(i));
  }
  
  arma::mat erTildeStar_mat(tmx, n, fill::zeros);
  arma::mat erTildeStarCent_mat(tmx, n, fill::zeros);
  
  for(int j=0; j<tmx; ++j){
    arma::vec betaTilde = trans(betaTilde_mat.row(j));
    
    arma::vec erTildeStar =
      trans(YTildeStar_mat.row(j)) - XTildeStar * betaTilde * scal;
    arma::vec erTildeStarCent = erTildeStar - mean(erTildeStar);
    erTildeStar_mat.row(j) = trans(erTildeStar);
    erTildeStarCent_mat.row(j) = trans(erTildeStarCent);
  }
  
  // centering
  arma::vec XTildeStarbar = trans(mean(XTildeStar));
  arma::vec YTildeStarbar_vec = mean(YTildeStar_mat, 1);
  arma::mat YTildeStarCent(tmx, n);
  
  arma::mat X0Cent_by_star(n0,tt); // center X0 by \bar{X}^*
  for(int j=0; j<tmx; ++j){
    YTildeStarCent.row(j) = YTildeStar_mat.row(j) - YTildeStarbar_vec(j);
  }
  arma::mat XTildeStarCent(n,tt);
  arma::mat DtTildeStar_mat(tmx, tt);
  for(int t=0;t<tt;++t){
    arma::vec XTildeStarCent_inner = XTildeStar.col(t) - XTildeStarbar(t);
    XTildeStarCent.col(t) = XTildeStarCent_inner;
    for(int j=0; j<tmx; ++j){
      DtTildeStar_mat(j, t) =
        mean(XTildeStarCent_inner %
        trans(YTildeStarCent.row(j)));
    }
    
    arma::vec X0Cent_by_star_inner = X0.col(t) - XTildeStarbar(t);
    X0Cent_by_star.col(t) = X0Cent_by_star_inner; // center X0 by \bar{X}^*
  }
  // X0 = X0Cent_by_star; // center X0 by \bar{X}^*
  
  
  // estimation
  
  arma::mat GaTildeStar = cov(XTildeStarCent); arma::vec gaTildeStar; arma::mat phiTildeStar;
  eig_sym(gaTildeStar, phiTildeStar, GaTildeStar);
  gaTildeStar = reverse(gaTildeStar) * scal;
  phiTildeStar = fliplr(phiTildeStar) / sqrt(scal);
  
  arma::cube GaInvTildeStar_cube(tt, tt, tmx, fill::zeros);
  arma::mat temppTildeStar(tt, tt, fill::zeros);
  for(int j=0; j<tmx; ++j){
    temppTildeStar = temppTildeStar +
      phiTildeStar.col(j) * trans(phiTildeStar.col(j)) / gaTildeStar(j);
    GaInvTildeStar_cube.slice(j) = temppTildeStar;
  }
  
  arma::cube betaTildeStar_cube_g(tmx, g_size, tt, fill::zeros);
  arma::cube Y0TildeStar_cube_g(tmx, g_size, n0, fill::zeros);
  arma::cube erTildeHatStar_cube_g(tmx, g_size, n, fill::zeros);
  arma::cube erTildeHatStarCent_cube_g(tmx, g_size, n, fill::zeros);
  arma::field<arma::mat> LdTildeStar_field_g(tmx, g_size);
  
  for(int ii_g=0; ii_g<g_size; ++ii_g){
    int g = gn_vec(ii_g);
    arma::vec UTilde = trans(UTilde_mat.row(g-1));
    arma::vec DtTildeStar = trans(DtTildeStar_mat.row(g-1));
    arma::vec YTildeStar = trans(YTildeStar_mat.row(g-1));
    
    for(int j=0; j<tmx; ++j){
      
      arma::vec betaTildeStar = GaInvTildeStar_cube.slice(j) *
        (DtTildeStar-UTilde) * scal;
      betaTildeStar_cube_g.tube(j, ii_g) = betaTildeStar;
      arma::vec Y0TildeStar_j = X0Cent_by_star * betaTildeStar * scal; // X0 appear
      Y0TildeStar_cube_g.tube(j, ii_g) = Y0TildeStar_j;
      
      arma::vec erTildeHatStar =
        YTildeStar - XTildeStar * betaTildeStar * scal;
      arma::vec erTildeHatStarCent = erTildeHatStar - mean(erTildeHatStar);
      erTildeHatStar_cube_g.tube(j, ii_g) = trans(erTildeHatStar);
      erTildeHatStarCent_cube_g.tube(j, ii_g) = trans(erTildeHatStarCent);
      
      arma::mat XTildeStarCent_erTildeHatStarCent(n, tt, fill::zeros);
      for(int i=0; i<n; ++i){
        XTildeStarCent_erTildeHatStarCent.row(i) =
          XTildeStarCent.row(i) * erTildeHatStarCent(i);
      }
      LdTildeStar_field_g(j, ii_g) = cov(XTildeStarCent_erTildeHatStarCent);
      
    }
  }
  
  // inference (testing only when enforcing the null)
  
  arma::field<arma::vec> statStar_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> statStar_std1_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> statStar_std2_field_inner(k_size, h_size, g_size);
  
  for(int ii_k=0; ii_k<k_size; ++ii_k){
    int k = kn_vec(ii_k);
    for(int ii_g=0; ii_g<g_size; ++ii_g){
      int g = gn_vec(ii_g);
      arma::mat LdTildeStar = LdTildeStar_field_g(k-1, ii_g);
      arma::vec Y0Tilde_g = trans(Y0Tilde_mat.row(g-1));
      
      for(int ii_h=0; ii_h<h_size; ++ii_h){
        int h = hn_vec(ii_h);
        
        arma::vec Y0TildeStar_h = Y0TildeStar_cube_g.tube(h-1, ii_g);
        arma::vec TTildeStar = Y0TildeStar_h - Y0Tilde_g;
        
        
        // std 1
        arma::vec sxHat_h = sxHat_cube_kh.tube(ii_k, ii_h) ; // X0 appear, but centered by \bar{X}
        arma::vec TTildeStar_std1 = sqrt(n) * TTildeStar / sqrt(sxHat_h);
        
        // std 2
        
        arma::mat GaInvTildeStar_h = GaInvTildeStar_cube.slice(h-1);
        arma::mat gixTildeStar = GaInvTildeStar_h * trans(X0Cent_by_star) * scal; // X0 appear
        arma::vec sxTildeStar_h =
          trans(sum((LdTildeStar * gixTildeStar * scal) % gixTildeStar * scal));
        
        arma::vec TTildeStar_std2 = sqrt(n) * TTildeStar / sqrt(sxTildeStar_h);
        
        
        
        // testing
        
        double statL2Star = sum(square(TTildeStar));
        double statMaxStar = max(abs(TTildeStar));
        statStar_field_inner(ii_k, ii_h, ii_g) =
          {statL2Star, statMaxStar};
        
        double statL2Star_std1 = sum(square(TTildeStar_std1));
        double statMaxStar_std1 = max(abs(TTildeStar_std1));
        statStar_std1_field_inner(ii_k, ii_h, ii_g) = 
          {statL2Star_std1, statMaxStar_std1};
        
        double statL2Star_std2 = sum(square(TTildeStar_std2));
        double statMaxStar_std2 = max(abs(TTildeStar_std2));
        statStar_std2_field_inner(ii_k, ii_h, ii_g) = 
          {statL2Star_std2, statMaxStar_std2};
        
        
      }
    }
  }
  
  List res_pb_inner = List::create(
    Named("statStar_field_inner") = statStar_field_inner,
    Named("statStar_std1_field_inner") = statStar_std1_field_inner,
    Named("statStar_std2_field_inner") = statStar_std2_field_inner
  );
  
  return(res_pb_inner);
}





//[[Rcpp::export]]
List PBouter(
    arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, 
    arma::mat& betaHat_mat, arma::mat& UHat_mat, arma::cube& sxHat_cube_kh,
    arma::mat& Y0Hat_mat, arma::cube& Y0Hat_trunc_hg_cube,
    arma::mat& YTilde_mat, arma::mat& betaTilde_mat, arma::mat& UTilde_mat, arma::mat& Y0Tilde_mat,
    arma::mat& stat_nonstd_mat, arma::cube& stat_Sstd_cube_kh,
    arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec,
    int M_bts, double ap = 0.05
){
  arma::vec ap_vec = {1-ap};
  arma::vec qp_vec = {ap*0.5, 1-ap*0.5};
  
  int n = X.n_rows; int n0 = X0.n_rows; int tt = X.n_cols;
  double scal = range(tGrid) / tt;
  int k_size = kn_vec.size(); 
  int h_size = hn_vec.size(); 
  int g_size = gn_vec.size(); 
  
  // arma::vec qp_vec = {ap*0.5, 1-ap*0.5};
  // arma::vec ap_vec = {1-ap};
  
  arma::vec tune(3); tune(0) = max(kn_vec); tune(1) = max(hn_vec); tune(2) = max(gn_vec);
  int tmx = tune.max();
  
  arma::field<arma::field<arma::vec>> rootsCI_trunc_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsCI_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsStd1CI_trunc_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsStd1CI_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsStd2CI_trunc_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsStd2CI_field(M_bts);
  
  arma::field<arma::field<arma::vec>> statStar_field(M_bts);
  arma::field<arma::field<arma::vec>> statStar_std1_field(M_bts);
  arma::field<arma::field<arma::vec>> statStar_std2_field(M_bts);
  arma::field<arma::field<arma::vec>> statTildeStar_field(M_bts);
  arma::field<arma::field<arma::vec>> statTildeStar_std1_field(M_bts);
  arma::field<arma::field<arma::vec>> statTildeStar_std2_field(M_bts);
  
  for(int m=0; m<M_bts; ++m){
    
    // resampling
    arma::vec idx = floor(Rcpp::runif(n)*n);
    
    List resPBinner1 = PBinner1(
      idx, X, Y, X0, tGrid, 
      betaHat_mat, UHat_mat, sxHat_cube_kh,
      Y0Hat_mat, Y0Hat_trunc_hg_cube,
      kn_vec, hn_vec, gn_vec
    );
    List resPBinner2 = PBinner2(
      idx, X, YTilde_mat, X0, tGrid, 
      betaTilde_mat, UTilde_mat, sxHat_cube_kh,
      Y0Tilde_mat, 
      kn_vec, hn_vec, gn_vec
    );
    
    arma::field<arma::vec> rootsCI_trunc_field_inner = 
      resPBinner1["rootsCI_trunc_field_inner"];
    arma::field<arma::vec> rootsCI_field_inner = 
      resPBinner1["rootsCI_field_inner"];
    arma::field<arma::vec> rootsStd1CI_trunc_field_inner = 
      resPBinner1["rootsStd1CI_trunc_field_inner"];
    arma::field<arma::vec> rootsStd1CI_field_inner = 
      resPBinner1["rootsStd1CI_field_inner"];
    arma::field<arma::vec> rootsStd2CI_trunc_field_inner = 
      resPBinner1["rootsStd2CI_trunc_field_inner"];
    arma::field<arma::vec> rootsStd2CI_field_inner = 
      resPBinner1["rootsStd2CI_field_inner"];
    
    rootsCI_trunc_field(m) = 
      rootsCI_trunc_field_inner;
    rootsCI_field(m) = 
      rootsCI_field_inner;
    rootsStd1CI_trunc_field(m) = 
      rootsStd1CI_trunc_field_inner;
    rootsStd1CI_field(m) = 
      rootsStd1CI_field_inner;
    rootsStd2CI_trunc_field(m) = 
      rootsStd2CI_trunc_field_inner;
    rootsStd2CI_field(m) = 
      rootsStd2CI_field_inner;
    
    arma::field<arma::vec> statStar_field_inner = 
      resPBinner1["statStar_field_inner"];
    arma::field<arma::vec> statStar_std1_field_inner = 
      resPBinner1["statStar_std1_field_inner"];
    arma::field<arma::vec> statStar_std2_field_inner = 
      resPBinner1["statStar_std2_field_inner"];
    
    arma::field<arma::vec> statTildeStar_field_inner = 
      resPBinner2["statStar_field_inner"];
    arma::field<arma::vec> statTildeStar_std1_field_inner = 
      resPBinner2["statStar_std1_field_inner"];
    arma::field<arma::vec> statTildeStar_std2_field_inner = 
      resPBinner2["statStar_std2_field_inner"];
    
    statStar_field(m) = 
      statStar_field_inner;
    statStar_std1_field(m) = 
      statStar_std1_field_inner;
    statStar_std2_field(m) = 
      statStar_std2_field_inner;
    
    statTildeStar_field(m) = 
      statTildeStar_field_inner;
    statTildeStar_std1_field(m) = 
      statTildeStar_std1_field_inner;
    statTildeStar_std2_field(m) = 
      statTildeStar_std2_field_inner;
    
  }
  
  
  arma::field<arma::cube> roots_field_khg(k_size, h_size, g_size);
  arma::field<arma::cube> intervals_field_khg(k_size, h_size, g_size);
  arma::field<arma::cube> stat_field_khg(k_size, h_size, g_size);
  arma::field<arma::mat> pVal_field_khg(k_size, h_size, g_size);
  
  for(int ii_k=0; ii_k<k_size; ++ii_k){
    int k = kn_vec(ii_k);
    for(int ii_h=0; ii_h<h_size; ++ii_h){
      int h = hn_vec(ii_h);
      arma::vec Y0Hat_h = trans(Y0Hat_mat.row(h-1));
      arma::vec sxHat = sxHat_cube_kh.tube(ii_k, ii_h);
      
      for(int ii_g=0; ii_g<g_size; ++ii_g){
        int g = gn_vec(ii_g);
        
        arma::cube roots_cube(6, M_bts, n0, fill::zeros);
        // Hat_nonstd, Hat_std1, Hat_std2,
        // Tilde_nonstd, Tilde_std1, Tilde_std2
        arma::cube stat_cube(6, M_bts, 2, fill::zeros);
        arma::cube pVal_cube(6, M_bts, 2, fill::zeros);
        
        for(int m=0; m<M_bts; ++m){
          
          // intervals without std
          
          arma::field<arma::vec> rootsCI_trunc_field_inner = rootsCI_trunc_field(m);
          arma::vec rootsCI_trunc_vec = rootsCI_trunc_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(0, m) = rootsCI_trunc_vec;
          
          arma::field<arma::vec> rootsCI_field_inner = rootsCI_field(m);
          arma::vec rootsCI_vec = rootsCI_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(1, m) = rootsCI_vec;
          
          // intervals with std1
          
          arma::field<arma::vec> rootsStd1CI_trunc_field_inner = rootsStd1CI_trunc_field(m);
          arma::vec rootsStd1CI_trunc_vec = rootsStd1CI_trunc_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(2, m) = rootsStd1CI_trunc_vec;
          
          arma::field<arma::vec> rootsStd1CI_field_inner = rootsStd1CI_field(m);
          arma::vec rootsStd1CI_vec = rootsStd1CI_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(3, m) = rootsStd1CI_vec;
          
          // intervals with std2
          
          arma::field<arma::vec> rootsStd2CI_trunc_field_inner = rootsStd2CI_trunc_field(m);
          arma::vec rootsStd2CI_trunc_vec = rootsStd2CI_trunc_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(4, m) = rootsStd2CI_trunc_vec;
          
          arma::field<arma::vec> rootsStd2CI_field_inner = rootsStd2CI_field(m);
          arma::vec rootsStd2CI_vec = rootsStd2CI_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(5, m) = rootsStd2CI_vec;
          
          // testing: not enforcing the null
          
          arma::field<arma::vec> statStar_field_inner = statStar_field(m);
          arma::vec statStar_vec = statStar_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(0, m) = statStar_vec;
          
          arma::field<arma::vec> statStar_std1_field_inner = statStar_std1_field(m);
          arma::vec statStar_std1_vec = statStar_std1_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(1, m) = statStar_std1_vec;
          
          arma::field<arma::vec> statStar_std2_field_inner = statStar_std2_field(m);
          arma::vec statStar_std2_vec = statStar_std2_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(2, m) = statStar_std2_vec;
          
          if(stat_nonstd_mat(ii_h, 0)<statStar_vec(0)){pVal_cube(0, m, 0)=1;}
          if(stat_nonstd_mat(ii_h, 1)<statStar_vec(1)){pVal_cube(0, m, 1)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 0)<statStar_std1_vec(0)){pVal_cube(1, m, 0)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 1)<statStar_std1_vec(1)){pVal_cube(1, m, 1)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 0)<statStar_std2_vec(0)){pVal_cube(2, m, 0)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 1)<statStar_std2_vec(1)){pVal_cube(2, m, 1)=1;}
          
          // testing: enforcing the null
          
          arma::field<arma::vec> statTildeStar_field_inner = statTildeStar_field(m);
          arma::vec statTildeStar_vec = statTildeStar_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(3, m) = statTildeStar_vec;
          
          arma::field<arma::vec> statTildeStar_std1_field_inner = statTildeStar_std1_field(m);
          arma::vec statTildeStar_std1_vec = statTildeStar_std1_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(4, m) = statTildeStar_std1_vec;
          
          arma::field<arma::vec> statTildeStar_std2_field_inner = statTildeStar_std2_field(m);
          arma::vec statTildeStar_std2_vec = statTildeStar_std2_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(5, m) = statTildeStar_std2_vec;
          
          if(stat_nonstd_mat(ii_h, 0)<statTildeStar_vec(0)){pVal_cube(3, m, 0)=1;}
          if(stat_nonstd_mat(ii_h, 1)<statTildeStar_vec(1)){pVal_cube(3, m, 1)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 0)<statTildeStar_std1_vec(0)){pVal_cube(4, m, 0)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 1)<statTildeStar_std1_vec(1)){pVal_cube(4, m, 1)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 0)<statTildeStar_std2_vec(0)){pVal_cube(5, m, 0)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 1)<statTildeStar_std2_vec(1)){pVal_cube(5, m, 1)=1;}
          
        }
        
        // ======
        // intervals
        
        roots_field_khg(ii_k, ii_h, ii_g) = roots_cube;
        arma::cube roots_cube_abs = abs(roots_cube);
        
        // individual intervals
        
        arma::cube II(n0, 2, 8);
        for(int l=0; l<2; ++l){
          // arma::mat temp = roots_cube.row(l);
          // arma::mat qrt = quantile(temp, qp_vec, 0);
          // arma::mat tempp = roots_cube_abs.row(l);
          // arma::mat qrt_sym = quantile(tempp, ap_vec, 0);
          for(int i0=0; i0<n0; ++i0){
            
            
            arma::mat temp = roots_cube.slice(i0);
            arma::rowvec tempp = temp.row(l);
            arma::rowvec qrt = quantile(tempp, qp_vec);
            
            arma::mat temp_abs = roots_cube_abs.slice(i0);
            arma::rowvec tempp_abs = temp_abs.row(l);
            arma::rowvec qrt_sym = quantile(tempp_abs, ap_vec);
            
            
            II(i0, 0, l) = Y0Hat_h(i0) - qrt(1);
            II(i0, 1, l) = Y0Hat_h(i0) - qrt(0);
            
            II(i0, 0, l+2) = Y0Hat_h(i0) - qrt_sym(0);
            II(i0, 1, l+2) = Y0Hat_h(i0) + qrt_sym(0);
            
          }
        }
        
        // II with std
        
        for(int l=4; l<6; ++l){
          // arma::mat temp = roots_cube.row(l);
          // arma::mat qrt = quantile(temp, qp_vec, 0);
          // arma::mat tempp = roots_cube_abs.row(l);
          // arma::mat qrt_sym = quantile(tempp, ap_vec, 0);
          for(int i0=0; i0<n0; ++i0){
            
            
            arma::mat temp = roots_cube.slice(i0);
            arma::rowvec tempp = temp.row(l);
            arma::rowvec qrt = quantile(tempp, qp_vec);
            
            arma::mat temp_abs = roots_cube_abs.slice(i0);
            arma::rowvec tempp_abs = temp_abs.row(l);
            arma::rowvec qrt_sym = quantile(tempp_abs, ap_vec);
            
            
            double seInv = sqrt( sxHat(i0) / n );
            
            II(i0, 0, l) = Y0Hat_h(i0) - seInv * qrt(1);
            II(i0, 1, l) = Y0Hat_h(i0) - seInv * qrt(0);
            
            II(i0, 0, l+2) = Y0Hat_h(i0) - seInv * qrt_sym(0);
            II(i0, 1, l+2) = Y0Hat_h(i0) + seInv * qrt_sym(0);
            
          }
        }
        
        // simultaneous intervals
        
        arma::vec qrt_max_rb(6, fill::zeros);
        for(int l=0; l<6; ++l){
          arma::vec temp = max(roots_cube_abs.row(l), 1);
          qrt_max_rb(l) = as_scalar(quantile(temp, ap_vec));
        }
        
        arma::cube SI(n0, 2, 6);
        for(int i0=0; i0<n0; ++i0){
          
          double seInv = sqrt( sxHat(i0) / n );
          
          for(int l=0; l<2; ++l){
            SI(i0, 0, l) = Y0Hat_h(i0) - qrt_max_rb(l);
            SI(i0, 1, l) = Y0Hat_h(i0) + qrt_max_rb(l);
            
            SI(i0, 0, l+2) = Y0Hat_h(i0) - qrt_max_rb(l+2) * seInv;
            SI(i0, 1, l+2) = Y0Hat_h(i0) + qrt_max_rb(l+2) * seInv;
            SI(i0, 0, l+4) = Y0Hat_h(i0) - qrt_max_rb(l+4) * seInv;
            SI(i0, 1, l+4) = Y0Hat_h(i0) + qrt_max_rb(l+4) * seInv;
            
          }
          
        }
        
        // store
        
        intervals_field_khg(ii_k, ii_h, ii_g) = join_slices(II, SI);
        
        // ======
        // testing
        
        stat_field_khg(ii_k, ii_h, ii_g) = stat_cube;
        
        pVal_field_khg(ii_k, ii_h, ii_g) = mean(pVal_cube, 1);
        
        
      }
    }
  }
  
  
  List resPBouter = List::create(
    Named("roots") = roots_field_khg,
    Named("intervals") = intervals_field_khg,
    Named("stat") = stat_field_khg,
    Named("pVal") = pVal_field_khg
  );
  
  return(resPBouter);
  
  
  
  
}

























//[[Rcpp::export]]
List WBinner(
    arma::vec& W, arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, 
    arma::mat& betaHat_mat, arma::mat& Y0Hat_mat, arma::cube& Y0Hat_trunc_hg_cube,
    arma::mat& erHat_mat, arma::mat& YHat_mat, arma::mat& XCent, 
    arma::mat& YTilde_mat, arma::mat& Y0Tilde_mat, arma::cube& GaInvHat_cube, arma::cube& sxHat_cube_kh,
    arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec
){
  
  int n = X.n_rows; int n0 = X0.n_rows; int tt = X.n_cols;
  double scal = range(tGrid) / tt;
  int k_size = kn_vec.size(); int h_size = hn_vec.size(); int g_size = gn_vec.size(); 
  arma::vec tune(3); tune(0) = max(kn_vec); tune(1) = max(hn_vec); tune(2) = max(gn_vec);
  int tmx = tune.max();
  
  
  // wild bootstrap
  
  arma::mat erStar_mat(k_size, n, fill::zeros);
  arma::mat erStarCent_mat(k_size, n, fill::zeros);
  
  arma::cube YStar_cube(k_size, g_size, n, fill::zeros);
  arma::cube YStarCent_cube(k_size, g_size, n, fill::zeros);
  
  arma::cube Y0Star_cube(k_size, g_size, n0, fill::zeros);
  
  // arma::field<arma::vec> Y0HatStar_field(k_size, h_size, g_size);
  
  arma::field<arma::vec> rootsCI_trunc_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsCI_field_inner(k_size, h_size, g_size);
  // arma::field<arma::vec> rootsPI_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsStdCI_trunc_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsStdCI_field_inner(k_size, h_size, g_size);
  // arma::field<arma::vec> rootsStdPI_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsStd2CI_trunc_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> rootsStd2CI_field_inner(k_size, h_size, g_size);
  
  arma::field<arma::vec> statStar_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> statStar_std_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> statStar_std2_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> statTildeStar_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> statTildeStar_std_field_inner(k_size, h_size, g_size);
  arma::field<arma::vec> statTildeStar_std2_field_inner(k_size, h_size, g_size);
  
  for(int ii_k=0; ii_k<k_size; ++ii_k){
    int k = kn_vec(ii_k);
    
    arma::vec erHat = trans(erHat_mat.row(k-1));
    arma::vec erStar(n);
    for(int i=0; i<n; ++i){ erStar(i) = erHat(i)*W(i); }
    arma::vec erStarCent = erStar - mean(erStar);
    
    erStar_mat.row(ii_k) = trans(erStar);
    erStarCent_mat.row(ii_k) = trans(erStarCent);
    
    // for LdStar
    
    arma::mat GaInvHat_k = GaInvHat_cube.slice(k-1);
    
    // // prediction
    // arma::vec idx0 = floor(Rcpp::runif(n0)*n);
    // arma::vec er0Star(n0);
    // for(int i0=0; i0<n0; ++i0){ er0Star(i0) = erHatCent(idx0(i0)); }
    
    for(int ii_g=0; ii_g<g_size; ++ii_g){
      int g = gn_vec(ii_g);
      
      arma::vec YHat_g = trans(YHat_mat.row(g-1));
      arma::vec YStar = YHat_g + erStar;
      double YStarbar = mean(YStar);
      arma::vec YStarCent = YStar - YStarbar;
      
      YStar_cube.tube(ii_k, ii_g) = YStar;
      YStarCent_cube.tube(ii_k, ii_g) = YStarCent;
      
      arma::vec DtHatStar(tt);
      for(int t=0;t<tt;++t){
        DtHatStar(t) = mean(XCent.col(t) % YStarCent);
      }
      
      
      
      // LdHatStar
      arma::vec betaHatStar_k = GaInvHat_k * DtHatStar * scal;
      vec erHatStar = YStar - X * betaHatStar_k * scal;
      vec erHatStarCent = erHatStar - mean(erHatStar);
      arma::mat XCent_erHatStarCent(n, tt, fill::zeros);
      for(int i=0; i<n; ++i){
        XCent_erHatStarCent.row(i) = XCent.row(i) * erHatStarCent(i);
      }
      arma::mat LdHatStar = cov(XCent_erHatStarCent);
      
      
      
      
      
      
      
      
      // enforcing the null
      
      arma::vec YTilde_g = trans(YTilde_mat.row(g-1));
      arma::vec YTildeStar = YTilde_g + erStar;
      double YTildeStarbar = mean(YTildeStar);
      arma::vec YTildeStarCent = YTildeStar - YTildeStarbar;
      
      YStar_cube.tube(ii_k, ii_g) = YStar;
      YStarCent_cube.tube(ii_k, ii_g) = YStarCent;
      
      arma::vec DtTildeStar(tt);
      for(int t=0;t<tt;++t){
        DtTildeStar(t) = mean(XCent.col(t) % YTildeStarCent);
      }
      
      arma::vec Y0Tilde_g = trans(Y0Tilde_mat.row(g-1));
      arma::vec Y0Hat_g = trans(Y0Hat_mat.row(g-1));
      
      // // prediction
      // arma::vec Y0Star = Y0Hat_g + er0Star;
      // Y0Star_cube.tube(ii_k, ii_g) = Y0Star;
      
      // arma::vec YHat_k = trans(YHat_mat.row(k));
      
      // LdTildeStar
      arma::vec betaTildeStar_k = GaInvHat_k * DtTildeStar * scal;
      arma::vec erTildeStar = YTildeStar - X * betaTildeStar_k * scal;
      arma::vec erTildeStarCent = erTildeStar - mean(erTildeStar);
      arma::mat XCent_erTildeStarCent(n, tt, fill::zeros);
      for(int i=0; i<n; ++i){
        XCent_erHatStarCent.row(i) = XCent.row(i) * erHatStarCent(i);
      }
      arma::mat LdTildeStar = cov(XCent_erHatStarCent);
      
      
      for(int ii_h=0; ii_h<h_size; ++ii_h){
        int h = hn_vec(ii_h);
        
        arma::mat GaInvHat_h = GaInvHat_cube.slice(h-1);
        arma::vec betaHatStar_h = GaInvHat_h * DtHatStar * scal;
        arma::vec Y0HatStar_h = X0 * betaHatStar_h * scal;
        
        // enforcing the null
        arma::vec betaTildeStar_h = GaInvHat_h * DtTildeStar * scal;
        arma::vec Y0TildeStar_h = X0 * betaTildeStar_h * scal;
        
        
        // Y0HatStar_field(ii_k, ii_h, ii_g) = Y0HatStar_field;
        
        // store for bootstrap intervals
        
        arma::vec Y0Hat_g_trunc_h = Y0Hat_trunc_hg_cube.tube(ii_h, ii_g);
        arma::vec rootsCI_trunc = Y0HatStar_h - Y0Hat_g_trunc_h;
        arma::vec TStar = Y0HatStar_h - Y0Hat_g;
        // arma::vec rootsPI = Y0HatStar_h - Y0Star;
        
        rootsCI_trunc_field_inner(ii_k, ii_h, ii_g) = rootsCI_trunc;
        rootsCI_field_inner(ii_k, ii_h, ii_g) = TStar;
        // rootsPI_field_inner(ii_k, ii_h, ii_g) = rootsPI;
        
        // arma::vec sxHat_h = trans(sxHat_cube_h.row(ii_h));
        
        
        
        
        
        
        
        
        
        
        
        // std 1
        
        arma::vec sxHat_h(n0);
        for(int i0=0; i0<n0; ++i0){
          sxHat_h(i0) = sxHat_cube_kh(ii_k,ii_h,i0);
        }
        
        rootsStdCI_trunc_field_inner(ii_k, ii_h, ii_g) =
          sqrt(n) * rootsCI_trunc / sqrt(sxHat_h);
        arma::vec TStar_std = sqrt(n) * TStar / sqrt(sxHat_h);
        rootsStdCI_field_inner(ii_k, ii_h, ii_g) = TStar_std;
        // rootsStdPI_field_inner(ii_k, ii_h, ii_g) =
        // sqrt(n) * rootsPI / sqrt(txHat_h/n+1);
        
        
        // std 2
        arma::mat gix = GaInvHat_h * trans(X0) * scal;
        arma::vec sxHatStar_h = trans(sum((LdHatStar * gix * scal) % gix * scal));
        
        rootsStd2CI_trunc_field_inner(ii_k, ii_h, ii_g) =
          sqrt(n) * rootsCI_trunc / sqrt(sxHatStar_h);
        arma::vec TStar_std2 = sqrt(n) * TStar / sqrt(sxHatStar_h);
        rootsStd2CI_field_inner(ii_k, ii_h, ii_g) = TStar_std2;
        
        
        
        // testing
        
        // not enforcing the null
        
        double statL2Star = sum(square(TStar));
        double statMaxStar = max(abs(TStar));
        statStar_field_inner(ii_k, ii_h, ii_g) = {statL2Star, statMaxStar};
        
        double statL2Star_std = sum(square(TStar_std));
        double statMaxStar_std = max(abs(TStar_std));
        statStar_std_field_inner(ii_k, ii_h, ii_g) = {statL2Star_std, statMaxStar_std};
        
        double statL2Star_std2 = sum(square(TStar_std2));
        double statMaxStar_std2 = max(abs(TStar_std2));
        statStar_std2_field_inner(ii_k, ii_h, ii_g) = {statL2Star_std2, statMaxStar_std2};
        
        // enforcing the null
        arma::vec TTildeStar = Y0TildeStar_h - Y0Tilde_g;
        double statTildeL2Star = sum(square(TTildeStar));
        double statTildeMaxStar = max(abs(TTildeStar));
        statTildeStar_field_inner(ii_k, ii_h, ii_g) =
          {statTildeL2Star, statTildeMaxStar};
        
        arma::vec TTildeStar_std = sqrt(n) * TTildeStar / sqrt(sxHat_h);
        double statTildeL2Star_std = sum(square(TTildeStar_std));
        double statTildeMaxStar_std = max(abs(TTildeStar_std));
        statTildeStar_std_field_inner(ii_k, ii_h, ii_g) =
          {statTildeL2Star_std, statTildeMaxStar_std};
        
        
        arma::vec sxTildeStar_h = trans(sum((LdTildeStar * gix * scal) % gix * scal));
        arma::vec TTildeStar_std2 = sqrt(n) * TTildeStar / sqrt(sxTildeStar_h);
        double statTildeL2Star_std2 = sum(square(TTildeStar_std2));
        double statTildeMaxStar_std2 = max(abs(TTildeStar_std2));
        statTildeStar_std2_field_inner(ii_k, ii_h, ii_g) =
          {statTildeL2Star_std2, statTildeMaxStar_std2};
        
        
        
      }
      
      
      
    }
  }
  
  List res_wb_inner = List::create(
    Named("rootsCI_trunc_field_inner") = rootsCI_trunc_field_inner,
    Named("rootsCI_field_inner") = rootsCI_field_inner,
    // Named("rootsPI_field_inner") = rootsPI_field_inner,
    Named("rootsStdCI_trunc_field_inner") = rootsStdCI_trunc_field_inner,
    Named("rootsStdCI_field_inner") = rootsStdCI_field_inner,
    // Named("rootsStdPI_field_inner") = rootsStdPI_field_inner,
    Named("rootsStd2CI_trunc_field_inner") = rootsStd2CI_trunc_field_inner,
    Named("rootsStd2CI_field_inner") = rootsStd2CI_field_inner,
    Named("statStar_field_inner") = statStar_field_inner,
    Named("statStar_std_field_inner") = statStar_std_field_inner,
    Named("statStar_std2_field_inner") = statStar_std2_field_inner,
    Named("statTildeStar_field_inner") = statTildeStar_field_inner,
    Named("statTildeStar_std_field_inner") = statTildeStar_std_field_inner,
    Named("statTildeStar_std2_field_inner") = statTildeStar_std_field_inner
  );
  
  // List res_wb_inner = List::create(X, Y);
  
  return(res_wb_inner);
  
}






//[[Rcpp::export]]
List WBouter(
    arma::mat& X, arma::vec& Y, arma::mat& X0, arma::vec& tGrid, 
    arma::mat& betaHat_mat, arma::mat& Y0Hat_mat, arma::cube& Y0Hat_trunc_hg_cube,
    arma::mat& erHat_mat, arma::mat& YHat_mat, arma::mat& XCent, 
    arma::mat& YTilde_mat, arma::mat& Y0Tilde_mat, arma::cube& GaInvHat_cube, arma::cube& sxHat_cube_kh,
    arma::mat& stat_nonstd_mat, arma::cube& stat_Sstd_cube_kh,
    arma::vec& kn_vec, arma::vec& hn_vec, arma::vec& gn_vec,
    arma::mat& W, double ap = 0.05){
  
  // W: n x M multiplier matrix
  
  arma::vec ap_vec = {1-ap};
  arma::vec qp_vec = {ap*0.5, 1-ap*0.5};
  
  int n = X.n_rows; int n0 = X0.n_rows; int tt = X.n_cols;
  double scal = range(tGrid) / tt;
  int k_size = kn_vec.size(); 
  int h_size = hn_vec.size(); 
  int g_size = gn_vec.size(); 
  
  // arma::vec qp_vec = {ap*0.5, 1-ap*0.5};
  // arma::vec ap_vec = {1-ap};
  
  arma::vec tune(3); tune(0) = max(kn_vec); tune(1) = max(hn_vec); tune(2) = max(gn_vec);
  int tmx = tune.max();
  int M_bts = W.n_cols;
  
  // 
  
  arma::field<arma::field<arma::vec>> rootsCI_trunc_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsCI_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsStdCI_trunc_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsStdCI_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsStd2CI_trunc_field(M_bts);
  arma::field<arma::field<arma::vec>> rootsStd2CI_field(M_bts);
  
  arma::field<arma::field<arma::vec>> statStar_field(M_bts);
  arma::field<arma::field<arma::vec>> statStar_std_field(M_bts);
  arma::field<arma::field<arma::vec>> statStar_std2_field(M_bts);
  arma::field<arma::field<arma::vec>> statTildeStar_field(M_bts);
  arma::field<arma::field<arma::vec>> statTildeStar_std_field(M_bts);
  arma::field<arma::field<arma::vec>> statTildeStar_std2_field(M_bts);
  
  
  
  
  
  
  
  
  
  // resampling
  
  for(int m=0; m<M_bts; ++m){
    
    arma::vec W_vec = W.col(m);
    
    List resWBinner = 
      WBinner(
        W_vec, X, Y, X0, tGrid, 
        betaHat_mat, Y0Hat_mat, Y0Hat_trunc_hg_cube,
        erHat_mat, YHat_mat, XCent, 
        YTilde_mat, Y0Tilde_mat, GaInvHat_cube, sxHat_cube_kh,
        kn_vec, hn_vec, gn_vec
      );
    
    
    
    arma::field<arma::vec> rootsCI_trunc_field_inner =
      resWBinner["rootsCI_trunc_field_inner"];
    arma::field<arma::vec> rootsCI_field_inner = 
      resWBinner["rootsCI_field_inner"];
    arma::field<arma::vec> rootsStdCI_trunc_field_inner = 
      resWBinner["rootsStdCI_trunc_field_inner"];
    arma::field<arma::vec> rootsStdCI_field_inner = 
      resWBinner["rootsStdCI_field_inner"];
    arma::field<arma::vec> rootsStd2CI_trunc_field_inner = 
      resWBinner["rootsStd2CI_trunc_field_inner"];
    arma::field<arma::vec> rootsStd2CI_field_inner = 
      resWBinner["rootsStd2CI_field_inner"];
    
    
    rootsCI_trunc_field(m) = rootsCI_trunc_field_inner;
    rootsCI_field(m) = rootsCI_field_inner;
    rootsStdCI_trunc_field(m) = rootsStdCI_trunc_field_inner;
    rootsStdCI_field(m) = rootsStdCI_field_inner;
    rootsStd2CI_trunc_field(m) = rootsStd2CI_trunc_field_inner;
    rootsStd2CI_field(m) = rootsStd2CI_field_inner;
    
    arma::field<arma::vec> statStar_field_inner = 
      resWBinner["statStar_field_inner"];
    arma::field<arma::vec> statStar_std_field_inner = 
      resWBinner["statStar_std_field_inner"];
    arma::field<arma::vec> statStar_std2_field_inner = 
      resWBinner["statStar_std2_field_inner"];
    arma::field<arma::vec> statTildeStar_field_inner = 
      resWBinner["statTildeStar_field_inner"];
    arma::field<arma::vec> statTildeStar_std_field_inner = 
      resWBinner["statTildeStar_std_field_inner"];
    arma::field<arma::vec> statTildeStar_std2_field_inner = 
      resWBinner["statTildeStar_std2_field_inner"];
    
    statStar_field(m) = statStar_field_inner;
    statStar_std_field(m) = statStar_std_field_inner;
    statStar_std2_field(m) = statStar_std2_field_inner;
    statTildeStar_field(m) = statTildeStar_field_inner;
    statTildeStar_std_field(m) = statTildeStar_std_field_inner;
    statTildeStar_std2_field(m) = statTildeStar_std2_field_inner;
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // inference
  
  arma::field<arma::cube> roots_field_khg(k_size, h_size, g_size);
  arma::field<arma::cube> intervals_field_khg(k_size, h_size, g_size);
  arma::field<arma::cube> stat_field_khg(k_size, h_size, g_size);
  arma::field<arma::mat> pVal_field_khg(k_size, h_size, g_size);
  
  for(int ii_k=0; ii_k<k_size; ++ii_k){
    int k = kn_vec(ii_k);
    for(int ii_h=0; ii_h<h_size; ++ii_h){
      int h = hn_vec(ii_h);
      arma::vec Y0Hat_h = trans(Y0Hat_mat.row(h-1));
      arma::vec sxHat_h = sxHat_cube_kh.tube(ii_k, ii_h) ;
      
      for(int ii_g=0; ii_g<g_size; ++ii_g){
        int g = gn_vec(ii_g);
        
        arma::cube roots_cube(6, M_bts, n0, fill::zeros);
        // Hat_nonstd, Hat_std, Tilde_nonstd, Tilde_std
        arma::cube stat_cube(6, M_bts, 2, fill::zeros);
        arma::cube pVal_cube(6, M_bts, 2, fill::zeros);
        
        
        for(int m=0; m<M_bts; ++m){
          
          // intervals without std
          
          arma::field<arma::vec> rootsCI_trunc_field_inner = rootsCI_trunc_field(m);
          arma::vec rootsCI_trunc_vec = rootsCI_trunc_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(0, m) = rootsCI_trunc_vec;
          
          arma::field<arma::vec> rootsCI_field_inner = rootsCI_field(m);
          arma::vec rootsCI_vec = rootsCI_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(1, m) = rootsCI_vec;
          
          // intervals with std
          
          arma::field<arma::vec> rootsStdCI_trunc_field_inner = rootsStdCI_trunc_field(m);
          arma::vec rootsStdCI_trunc_vec = rootsStdCI_trunc_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(2, m) = rootsStdCI_trunc_vec;
          
          arma::field<arma::vec> rootsStdCI_field_inner = rootsStdCI_field(m);
          arma::vec rootsStdCI_vec = rootsStdCI_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(3, m) = rootsStdCI_vec;
          
          // intervals with std2
          
          arma::field<arma::vec> rootsStd2CI_trunc_field_inner = rootsStd2CI_trunc_field(m);
          arma::vec rootsStd2CI_trunc_vec = rootsStd2CI_trunc_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(4, m) = rootsStd2CI_trunc_vec;
          
          arma::field<arma::vec> rootsStd2CI_field_inner = rootsStd2CI_field(m);
          arma::vec rootsStd2CI_vec = rootsStd2CI_field_inner(ii_k, ii_h, ii_g);
          roots_cube.tube(5, m) = rootsStd2CI_vec;
          
          // testing: not enforcing the null
          
          arma::field<arma::vec> statStar_field_inner = statStar_field(m);
          arma::vec statStar_vec = statStar_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(0, m) = statStar_vec;
          
          arma::field<arma::vec> statStar_std_field_inner = statStar_std_field(m);
          arma::vec statStar_std_vec = statStar_std_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(1, m) = statStar_std_vec;
          
          arma::field<arma::vec> statStar_std2_field_inner = statStar_std2_field(m);
          arma::vec statStar_std2_vec = statStar_std2_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(2, m) = statStar_std2_vec;
          
          if(stat_nonstd_mat(ii_h, 0)<statStar_vec(0)){pVal_cube(0, m, 0)=1;}
          if(stat_nonstd_mat(ii_h, 1)<statStar_vec(1)){pVal_cube(0, m, 1)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 0)<statStar_std_vec(0)){pVal_cube(1, m, 0)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 1)<statStar_std_vec(1)){pVal_cube(1, m, 1)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 0)<statStar_std2_vec(0)){pVal_cube(2, m, 0)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 1)<statStar_std2_vec(1)){pVal_cube(2, m, 1)=1;}
          
          // testing: enforcing the null
          
          arma::field<arma::vec> statTildeStar_field_inner = statTildeStar_field(m);
          arma::vec statTildeStar_vec = statTildeStar_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(3, m) = statTildeStar_vec;
          
          arma::field<arma::vec> statTildeStar_std_field_inner = statTildeStar_std_field(m);
          arma::vec statTildeStar_std_vec = statTildeStar_std_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(4, m) = statTildeStar_std_vec;
          
          arma::field<arma::vec> statTildeStar_std2_field_inner = statTildeStar_std2_field(m);
          arma::vec statTildeStar_std2_vec = statTildeStar_std2_field_inner(ii_k, ii_h, ii_g);
          stat_cube.tube(5, m) = statTildeStar_std2_vec;
          
          if(stat_nonstd_mat(ii_h, 0)<statTildeStar_vec(0)){pVal_cube(3, m, 0)=1;}
          if(stat_nonstd_mat(ii_h, 1)<statTildeStar_vec(1)){pVal_cube(3, m, 1)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 0)<statTildeStar_std_vec(0)){pVal_cube(4, m, 0)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 1)<statTildeStar_std_vec(1)){pVal_cube(4, m, 1)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 0)<statTildeStar_std2_vec(0)){pVal_cube(5, m, 0)=1;}
          if(stat_Sstd_cube_kh(ii_k, ii_h, 1)<statTildeStar_std2_vec(1)){pVal_cube(5, m, 1)=1;}
        }
        
        // ======
        // intervals
        
        roots_field_khg(ii_k, ii_h, ii_g) = roots_cube;
        arma::cube roots_cube_abs = abs(roots_cube);
        
        // individual intervals
        
        arma::cube II(n0, 2, 8);
        for(int l=0; l<2; ++l){
          // arma::mat temp = roots_cube.row(l);
          // arma::mat qrt = quantile(temp, qp_vec, 0);
          // arma::mat tempp = roots_cube_abs.row(l);
          // arma::mat qrt_sym = quantile(tempp, ap_vec, 0);
          for(int i0=0; i0<n0; ++i0){
            
            
            arma::mat temp = roots_cube.slice(i0);
            arma::rowvec tempp = temp.row(l);
            arma::rowvec qrt = quantile(tempp, qp_vec);
            
            arma::mat temp_abs = roots_cube_abs.slice(i0);
            arma::rowvec tempp_abs = temp_abs.row(l);
            arma::rowvec qrt_sym = quantile(tempp_abs, ap_vec);
            
            
            II(i0, 0, l) = Y0Hat_h(i0) - qrt(1);
            II(i0, 1, l) = Y0Hat_h(i0) - qrt(0);
            
            II(i0, 0, l+2) = Y0Hat_h(i0) - qrt_sym(0);
            II(i0, 1, l+2) = Y0Hat_h(i0) + qrt_sym(0);
            
          }
        }
        
        // II with std
        
        for(int l=4; l<6; ++l){
          // arma::mat temp = roots_cube.row(l);
          // arma::mat qrt = quantile(temp, qp_vec, 0);
          // arma::mat tempp = roots_cube_abs.row(l);
          // arma::mat qrt_sym = quantile(tempp, ap_vec, 0);
          for(int i0=0; i0<n0; ++i0){
            
            
            arma::mat temp = roots_cube.slice(i0);
            arma::rowvec tempp = temp.row(l);
            arma::rowvec qrt = quantile(tempp, qp_vec);
            
            arma::mat temp_abs = roots_cube_abs.slice(i0);
            arma::rowvec tempp_abs = temp_abs.row(l);
            arma::rowvec qrt_sym = quantile(tempp_abs, ap_vec);
            
            
            double seInv = sqrt( sxHat_h(i0) / n );
            
            II(i0, 0, l) = Y0Hat_h(i0) - seInv * qrt(1);
            II(i0, 1, l) = Y0Hat_h(i0) - seInv * qrt(0);
            
            II(i0, 0, l+2) = Y0Hat_h(i0) - seInv * qrt_sym(0);
            II(i0, 1, l+2) = Y0Hat_h(i0) + seInv * qrt_sym(0);
            
          }
        }
        
        // simultaneous intervals
        
        arma::vec qrt_max_rb(6, fill::zeros);
        for(int l=0; l<6; ++l){
          arma::vec temp = max(roots_cube_abs.row(l), 1);
          qrt_max_rb(l) = as_scalar(quantile(temp, ap_vec));
        }
        
        arma::cube SI(n0, 2, 6);
        for(int i0=0; i0<n0; ++i0){
          
          double seInv = sqrt( sxHat_h(i0) / n );
          
          for(int l=0; l<2; ++l){
            SI(i0, 0, l) = Y0Hat_h(i0) - qrt_max_rb(l);
            SI(i0, 1, l) = Y0Hat_h(i0) + qrt_max_rb(l);
            
            SI(i0, 0, l+2) = Y0Hat_h(i0) - qrt_max_rb(l+2) * seInv;
            SI(i0, 1, l+2) = Y0Hat_h(i0) + qrt_max_rb(l+2) * seInv;
            SI(i0, 0, l+4) = Y0Hat_h(i0) - qrt_max_rb(l+4) * seInv;
            SI(i0, 1, l+4) = Y0Hat_h(i0) + qrt_max_rb(l+4) * seInv;
            
          }
          
        }
        
        // store
        
        intervals_field_khg(ii_k, ii_h, ii_g) = join_slices(II, SI);
        
        // ======
        // testing
        
        stat_field_khg(ii_k, ii_h, ii_g) = stat_cube;
        
        pVal_field_khg(ii_k, ii_h, ii_g) = mean(pVal_cube, 1);
        
        
      }
    }
  }
  
  List resWBouter = List::create(
    Named("roots") = roots_field_khg,
    Named("intervals") = intervals_field_khg,
    Named("stat") = stat_field_khg,
    Named("pVal") = pVal_field_khg
  );
  
  // List resWBouter = List::create(rootsCI_trunc_field);
  
  return(resWBouter);
  
}
















