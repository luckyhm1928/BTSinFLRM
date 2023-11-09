#' @title Multiplier wild bootstrap in functional linear regression models
#' 
#' @description 
#' Conduct wild bootstrap inference for projections of new regressors onto the slope function 
#' in functional linear regression models with scalar response under heteroscedasticity.
#'
#' Apply the central limit theorem (CLT) (Yeon, Dai, and Nordman, 2024a) to find individual confidence intervals.
#' Apply the multiplier wild bootstrap (Yeon, Dai, and Nordman, 2024b) to find individual and simultaneous intervals.
#' Intervals are provided with or without either symmetrization or studentization.
#' Simultaneous intervals are studentized by either the estimated scaling or the bootstrap scaling.
#' Apply the bootstrap hypothesis tests 
#' of whether the projections onto multiple new regressors are simultaneously zero or not;
#' this testing procedure is similar to the ones proposed by (Yeon, Dai, and Nordman, 2024a).
#'
#' @param X An n by p matrix of regressor curves. Each row represents one observed regressor curve.
#' @param Y An vector of n responses. 
#' @param X0 An n0 by p matrix of new regressor curves. Each row represents one new regressor curve.
#' @param tGrid A vector of p densely equi-spaced time grid points where the curve are evaluated.
#' @param kn_vec A vector of truncation level k associated with the estimated error variance.
#' @param hn_vec A vector of truncation level h associated with estimation of the slope function.
#' @param gn_vec A vector of truncation level g where the estimated slope function plays a role of the true parameter in the bootstrap world.
#' @param M_bts An integer. The number of bootstrap resamples with default 1000.
#' @param multiplier A string. The type of multipliers used in the bootstrap procedure. 
#' Only one of the following three types is allowed: bin, normal, and HM. See the details.
#' @param ap A numeric scalar. The significance level for inference.
#' @return A list containing the following fields:
#' \item{est}{The FPCR estimated slope function with different trunction levels} 
#' \item{CLT}{The individual intervals based on CLT with different trunction levels}   
#' \item{WBintervals}{The intervals based on wild bootstrap with different truncation levels. 
#' Intervals are provided with or without either symmetrization or studentization.}
#' \item{WBpvalues}{P-values from the bootstrap hypothesis tests 
#' of whether the projections onto multiple new regressors are simultaneously zero or not.}
#' 
#' @details 
#' 
#' The following three types of multipliers are considered: 
#' two-point distribution with \eqn{\mathsf{P}(W_i = -(\sqrt{5}-1)/2) = (\sqrt{5}+1) / (2\sqrt{5}) = 1-\mathsf{P}(W_i = (\sqrt{5}+1)/2)}, 
#' the standard normal distribution with \eqn{W_i \sim \mathsf{N}(0,1)}, 
#' and the distribution of \eqn{W_i = V_i/2+(V_i^2-1)/2} 
#' for independent standard normal variables \eqn{\{V_i\}_{i=1}^n}, 
#' which gives the third moment one \eqn{\mathsf{E}^*[W_i^3]=1}.
#' These respectively correspond to \code{multiplier="bin"}, \code{multiplier="normal"}, and \code{multiplier="HM"}.
#' 
#' @seealso
#' \code{\link{RBinFLRM}}
#' \code{\link{PBinFLRM}}
#' 
#' @examples
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
#' # wild bootstrap with normal multipliers
#' M_bts=100
#' WBinFLRM(X, Y, X0, tGrid, kn_vec, hn_vec, gn_vec, M_bts, multiplier = "normal")
#' 
#' @export
WBinFLRM = function(X, Y, X0, tGrid, kn_vec, hn_vec, gn_vec, M_bts=1000, ap=0.05, multiplier="normal"){
  n0 = nrow(X0)
  res = WBinFLRM::inferFLRM(X, Y, X0, tGrid, kn_vec, hn_vec, gn_vec, ap)
  
  # CLT
  CLT_hetero = abind::abind(
    res$CLThetero$lowerCI_CLT_hetero,
    res$CLThetero$upperCI_CLT_hetero,
    along=4)
  dimnames(CLT_hetero) = list(
    paste0("kn=", kn_vec),
    paste0("hn=", hn_vec),
    paste0("X0_", 1:n0),
    c("lower", "upper")
  )
  CLT_hetero = aperm(CLT_hetero, c(2,4,3,1))
  
  # PB
  
  n=length(Y)
  p0 = (sqrt(5)+1) / (2*sqrt(5))
  
  if(multiplier == "normal"){
    Worigin = rnorm(n*M_bts)
  }else if(multiplier == "bin"){
    Worigin = (2*rbinom(n*M_bts, 1, p0)-1)*sqrt(5)/2+1/2
  }else if(multiplier == "HM"){
    V = rnorm(n*M_bts)
    Worigin = V/2 + (V^2-1) / 2
  }else{
    stop("Multiplier should be one of three types: bin, normal, and HM.")
  }
  W = matrix(Worigin, nrow=n)
  
  XCent = t(t(X) - colMeans(X))
  resWBouter = WBinFLRM::WBouter(
    X, Y, X0, tGrid, 
    res$est$betaHat, res$est$Y0Hat, res$est$Y0Hat_trunc_hg,
    res$est$erHat, res$est$YHat, XCent, 
    res$estTilde$YTilde, res$estTilde$Y0Tilde,
    res$est$GaInvHat, res$CLThetero$sxHat, 
    res$STAThetero$stat_nonstd, res$STAThetero$stat_Sstd,
    kn_vec, hn_vec, gn_vec, W, ap
    )
  dimnames(resWBouter$intervals) = list(
    paste0("kn=", kn_vec),
    paste0("hn=", hn_vec),
    paste0("gn=", gn_vec)
  )
  dimnames(resWBouter$pVal) = list(
    paste0("kn=", kn_vec),
    paste0("hn=", hn_vec),
    paste0("gn=", gn_vec)
  )
  
  for(ii_k in 1:length(kn_vec)){
    for(ii_h in 1:length(hn_vec)){
      for(ii_g in 1:length(gn_vec)){
        # intervals
        ress = abind::abind(
          array( resWBouter$intervals[[ii_k, ii_h, ii_g]][,,1:2], dim = c(n0, 2, 2) ),
          array( resWBouter$intervals[[ii_k, ii_h, ii_g]][,,3:4], dim = c(n0, 2, 2) ),
          array( resWBouter$intervals[[ii_k, ii_h, ii_g]][,,5:6], dim = c(n0, 2, 2) ),
          array( resWBouter$intervals[[ii_k, ii_h, ii_g]][,,7:8], dim = c(n0, 2, 2) ),
          array( resWBouter$intervals[[ii_k, ii_h, ii_g]][,,9:10], dim = c(n0, 2, 2) ),
          array( resWBouter$intervals[[ii_k, ii_h, ii_g]][,,11:12], dim = c(n0, 2, 2) ),
          array( resWBouter$intervals[[ii_k, ii_h, ii_g]][,,13:14], dim = c(n0, 2, 2) ),
          along=4
        )
        dimnames(ress) = list(
          paste0("X0_", 1:n0),
          c("lower", "upper"),
          c("trunc", "untrunc"),
          c("II_unsym/nonstd", "II_sym/nonstd", 
            "II_unsym/std",    "II_sym/std",  
            "SI_nonstd", "SI_std1", "SI_std2")
        )
        
        resWBouter$intervals[[ii_k,ii_h,ii_g]] = ress
        
        
        # p-values
        resss = abind::abind(
          resWBouter$pVal[[ii_k,ii_h,ii_g]][1:3,],
          resWBouter$pVal[[ii_k,ii_h,ii_g]][4:6,],
          along=3
        )
        dimnames(resss) = list(
          c("nonstd", "std1", "std2"),
          c('L2', "max"),
          c("hat", "tilde")
        )
        resWBouter$pVal[[ii_k,ii_h,ii_g]] = resss
        
      }
    }
  }
  
  # return
  resWB = list(
    est = res$est$betaHat,
    CLT = CLT_hetero,
    WBintervals = resWBouter$intervals,
    WBpvalues = resWBouter$pVal
  )
  return(resWB)
  
}
