#' @title Residual bootstrap in functional linear regression models
#' 
#' @description 
#' Find residual bootstrap intervals for projections of new regressors onto the slope function 
#' in functional linear regression models with scalar response under homoscedasticity.
#'
#' Apply the central limit theorem (CLT) (Yeon, Dai, and Nordman, 2023) to find individual confidence and prediction intervals.
#' Apply the residual bootstrap (Yeon, Dai, and Nordman, 2023) to find individual and simultaneous confidence/prediction intervals.
#' Individual intervals are provided with or without symmetrization.
#' Simultaneous intervals are provided with or without studentization.
#'
#' @param X An n by p matrix of regressor curves. Each row represents one observed regressor curve.
#' @param Y An vector of n responses. 
#' @param X0 An n0 by p matrix of new regressor curves. Each row represents one new regressor curve.
#' @param tGrid A vector of p densely equi-spaced time grid points where the curve are evaluated.
#' @param kn_vec A vector of truncation level k associated with the estimated error variance.
#' @param hn_vec A vector of truncation level h associated with estimation of the slope function.
#' @param gn_vec A vector of truncation level g where the estimated slope function plays a role of the true parameter in the bootstrap world.
#' @param M_bts An integer. The number of bootstrap resamples with default 1000.
#' @param ap A numeric scalar. The significance level for inference.
#' @return A list containing the following fields:
#' \item{est}{The FPCR estimated slope function with different trunction levels} 
#' \item{CLT}{The individual intervals based on CLT with different trunction levels}   
#' \item{RBintervals}{The intervals based on residual bootstrap with different truncation levels. 
#' Individual intervals are provided with or without symmetrization. 
#' Simultaneous intervals are provided with or without studentization by the estimated scaling.}
#' 
#' @seealso
#' \code{\link{PBinFLRM}}
#' \code{\link{WBinFLRM}}
#' 
#' @examples
#' 
#' # example
#' 
#' library(BTSinFLRM)
#' 
#' # simulate the regressor curves from Karhunen-Loeve expansion and responses
#' J=15
#' 
#' # eigenvalues
#' c=2; a=3.5; dt = c*(1:J)^(-a)
#' g1 = c*VGAM::zeta(a)
#' ga = c(g1, sapply(1:(J-1), function(j){g1 - sum(dt[1:j])}))
#' 
#' # eigenfunctions
#' tt = 100; tGrid = seq(0, 1, len=tt); J=15
#' phi= t(fda::fourier(tGrid, J))
#' 
#' # slope function
#' b=3.5; bj = c*((1:J)^(-b))*(rbinom(J, 1, 0.5)*2-1)
#' beta = c(bj %*% phi)
#' 
#' # generate (Gaussian) regressors and responses
#' n=50; n0=3
#' xi = matrix(rnorm(n*J), n, J); X = xi %*% (phi * sqrt(ga))
#' xi0 = matrix(rnorm(n0*J), n0, J); X0 = xi0 %*% (phi * sqrt(ga))
#' er = runif(n); Y = c(X %*% beta / tt + er)
#' 
#' # trunction parameters
#' kn_vec = c(3,4); hn_vec = c(4,5,6); gn_vec = c(5)
#' 
#' # residual bootstrap
#' M_bts=100
#' RBinFLRM(X, Y, X0, tGrid, kn_vec, hn_vec, gn_vec, M_bts)
#' 
#' @export
RBinFLRM = function(X, Y, X0, tGrid, kn_vec, hn_vec, gn_vec, M_bts=1000, ap=0.05){
  
  # center X0 by \bar{X}
  Xmean = colMeans(X)
  X0Cent = t(t(X0) - Xmean)
  
  # inference
  res = BTSinFLRM:::inferFLRM(X, Y, X0Cent, tGrid, kn_vec, hn_vec, gn_vec, ap)
  n0=nrow(X0)
  # CLT
  CI_CLT_homo = abind::abind(
    res$CLThomo$lowerCI_CLT_homo,
    res$CLThomo$upperCI_CLT_homo,
    along=4)
  PI_CLT_homo = abind::abind(
    res$CLThomo$lowerPI_CLT_homo,
    res$CLThomo$upperPI_CLT_homo,
    along=4)
  CLT_homo = abind::abind(CI_CLT_homo, PI_CLT_homo, along=5)
  dimnames(CLT_homo) = list(
    paste0("kn=", kn_vec),
    paste0("hn=", hn_vec),
    paste0("X0_", 1:n0),
    c("lower", "upper"),
    c("CI", "PI")
  )
  CLT_homo = aperm(CLT_homo, c(2,4,3,5,1))
  
  # RB
  XCent = t(t(X) - Xmean)
  resRBouter = BTSinFLRM:::RBouter(
    X, Y, X0Cent, tGrid, 
    res$est$betaHat, res$est$Y0Hat, res$est$Y0Hat_trunc_hg,
    res$est$erHatCent, res$est$YHat, XCent, 
    res$estTilde$YTilde, res$estTilde$Y0Tilde,
    res$est$GaInvHat, res$CLThomo$txHat, 
    res$STAThomo$stat_nonstd, res$STAThomo$stat_Tstd,
    kn_vec, hn_vec, gn_vec, M_bts, ap)
  dimnames(resRBouter$intervals) = list(
    paste0("kn=", kn_vec),
    paste0("hn=", hn_vec),
    paste0("gn=", gn_vec)
  )
  namesInt = c("CI_trunc", "CI", "PI")
  for(ii_k in 1:length(kn_vec)){
    for(ii_h in 1:length(hn_vec)){
      for(ii_g in 1:length(gn_vec)){
        dimnames(resRBouter$intervals[[ii_k,ii_h,ii_g]]) = list(
          paste0("X0_", 1:n0),
          c("lower", "upper"),
          c(
            paste0("I", namesInt, "_unsym"),
            paste0("I", namesInt, "_sym"),
            paste0("S", namesInt, "_nonstd"),
            paste0("S", namesInt, "_std")
          )
        )
      }
    }
  }
  
  # return
  resRB = list(
    est = res$est$betaHat,
    CLT = CLT_homo,
    RBintervals = resRBouter$intervals
  )
  
  return(resRB)
}
