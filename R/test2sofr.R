#' @title Bootstrap tests for two-sample functional linear regression
#' 
#' @description 
#' Conduct bootstrap tests by Yeon and Kokoszka (2025) for assessing the equality of the slope functions in two scalar-on-function linear regression models.
#' The test statistics are built based on either eigendecomposition of the pooled covariance operator (PD) or simultaneously diagonalized components (SD) from group-specific covariance operators.
#' Residual-type bootstrap methods are employed to improve accuracy either imposing the null or not.
#' The truncation levels are chosen by using fraction of variance explained (FVE) or relative variance explained (RVE).
#'
#' @param Y A vector of n responses from two gruops. 
#' @param X An n by p matrix of regressor curves from two groups. Each row represents one observed regressor curve.
#' @param tGrid A vector of p densely equi-spaced time grid points where the curve are evaluated.
#' @param group A vector of n elements that indicates two groups. 
#' @param h_max A positive integer value that represents the maximum of candidate truncation levels. 
#' @param rho_vec A vector of the threshold values using to choose the truncation levels; these must be numeric values between 0 and 1.
#' @param B A positive integer. The number of bootstrap resamples with default 1000.
#' @return A list containing the following fields:
#' \item{res_hEach}{The test statistics, the corresponding p-values, and the values of variance explained for each truncation level from 1 to h_max} 
#' \item{res_hSel}{The test statistics, the corresponding p-values, and the selected truncation levels by using different types of variance explained with threshold values in rho_vec}   
#' 
#' @references 
#' 
#' Yeon, H. and Kokoszka, P. (2025). Bootstrap two-sample tests for scalar-on-function regression. \emph{Submitted}
#' 
#' @examples
#' 
#' # example
#' 
#' library(BTSinFLRM)
#' 
#' Jtrue=20 # The number of basis elements used in data generation
#'
#' # create eigenvalues based on polynomial eigengaps with rate a
#' ev_poly = function(a, c, J){
#'   dt = c*(1:J)^(-a)
#'   ld1 = c*VGAM::zeta(a)
#'   ld = c(ld1, sapply(1:(J-1), function(j){ld1 - sum(dt[1:j])}))
#'   return(ld)
#' }
#'
#' # eigenvalues for unequal cov
#' ga_rough  = ev_poly(2.5, 2, Jtrue) # a=2.5
#' ga_smooth = ev_poly(5, 2, Jtrue)   # a=5
#'
#' # Equi-spaced time grid points in [0,1]
#' tt = 50; tGrid_ttt = seq(0, 1, len=tt+1); 
#' tGrid = tGrid_ttt[1:tt]
#' diffrange = diff(tGrid)[1] + diff(range(tGrid)) # when using left points
#' scal = diffrange / tt   # scaling factor for integration
#'
#' # basis functions for beta, and for X under equal cov
#' phi_base = t(fda::fourier(tGrid_ttt, 2*Jtrue))[-1,1:tt]  # trigonometric functions
#' phi = phi_base[1:Jtrue,]
#'
#' # basis functions for X under unequal cov
#' mono = t(sapply(1:Jtrue, function(j)tGrid^j))
#' phi_mono = BTSinFLRM:::orthoL2Equidense(mono, tGrid)
#' cheb = t(sapply(1:Jtrue, function(j)pracma::chebPoly(j, 2*tGrid-1)))
#' phi_cheb = BTSinFLRM:::orthoL2Equidense(cheb, tGrid)
#'
#' # Slope functions
#' b=2; Wbetaj = rbinom(Jtrue, 1, 0.5)*2-1
#' bj = 3*((1:Jtrue)^(-b))*Wbetaj; beta = c(bj %*% phi)
#'
#' beta_arr = array(0, dim=c(2, tt, 2), dimnames=list(
#'   1:2, NULL, paste0("type", 0:1)
#' ))
#' # 0: linear difference
#' beta_arr[1,,"type0"] = 1*beta
#' beta_arr[2,,"type0"] = 2*beta
#'
#' # 1: orthogonal difference
#' beta_arr[1,,"type1"] = c(bj[1:5] %*% phi[1:5,])
#' beta_arr[2,,"type1"] = c(bj[1:5] %*% phi[6:10,])
#'
#' beta1 = beta_arr[1,,"type0"]
#' c_alter = 0.5 # strength of alternative
#' beta2 = (1-c_alter)*beta1 + c_alter*beta_arr[2,,"type0"]
#'
#' # Determine the second order structures
#' eveq = TRUE # eigenvalues are equal?
#' efeq = TRUE # eigenfunctions are equal?
#' if(eveq){  # equal eigenvalues
#'   ga1 = ga2 = ga_rough
#' }else{  # unequal eigenvalues
#'   ga1 = ga_smooth; ga2 = ga_rough
#' }
#' if(efeq){  # equal eigenfunctions
#'   phi1 = phi; phi2 = phi
#' }else{  # unequal eigenfunctions
#'   phi1 = phi_mono; phi2 = phi_cheb
#' }
#'
#' # Other factors
#'
#' n=50; n1 = n2 = n/2      # n1=48; n2=52
#' group = c(rep(1,n1), rep(0,n2))
#' idx1 = as.logical(group); idx2 = !idx1
#' xi_type=0 # Gaussian random function
#' er_type=0 # Normal error
#' sig_sq =1 # the error variance
#'
#' # # # # # # # # # # # 
#' # generate data
#'
#' # FPC scores
#' Wj = matrix(rnorm(n*Jtrue), n, Jtrue)
#' if(xi_type==2){
#'   xi = rexp(n) - 1 # normal*exp FPC scores
#' }else if(xi_type==1){
#'   xi = rnorm(n) # normal*normal FPC scores
#' }else{
#'   xi = rep(1,n) # independent FPC scores
#' }
#'
#' # # unequal FPC scores
#' # xi1 = rnorm(n1); xi2 = rexp(n2) - 1
#' # xi[idx1,] = xi1; xi[idx2,] = xi2
#'
#' Xfpc = matrix(0, n, Jtrue)
#' Xfpc[idx1,] = t(sqrt(ga1) * t(xi[idx1]*Wj[idx1,]))
#' Xfpc[idx2,] = t(sqrt(ga2) * t(xi[idx2]*Wj[idx2,]))
#'
#' # regerssor functions
#'
#' X = matrix(0, n, tt)
#' X[idx1,] = Xfpc[idx1,] %*% phi1
#' X[idx2,] = Xfpc[idx2,] %*% phi2
#'
#' X1 = X[idx1,]; X2 = X[idx2,]
#'
#' X1beta1 = c(X1 %*% beta1 * scal)
#' X2beta2 = c(X2 %*% beta2 * scal)
#'
#' if(er_type==1){   # centered exponential error
#'   par_scale = sqrt(sig_sq)
#'   er = rexp(n, 1/par_scale) - par_scale
#' }else{     # normal error
#'   er = rnorm(n, 0,sd=sqrt(sig_sq))          
#' }
#'
#' # response variables
#' Y1 = X1beta1 + er[idx1]
#' Y2 = X2beta2 + er[idx2]
#' Y = rep(0,n); Y[idx1] = Y1; Y[idx2] = Y2
#'
#' # Bootstrap tests for two-sample scalar-on-function regression (not run)
#' # res = test2sofr(Y, X, tGrid, group)
#' 
#' @export
test2sofr = function(
    Y, X, tGrid, group, 
    h_max=20, 
    rho_vec = c(seq(0.75, 0.95, by=0.05), 0.99, 0.995, 0.999),
    B=1000, ap=0.05
    ){
  
  
  # =================est================
  # estimation and test statsitc
  
  
  n=nrow(X); tt=ncol(X); # tGrid=1:tt/tt
  diffrange = diff(tGrid)[1] + diff(range(tGrid)) # when using left points
  scal = diffrange / tt   # scaling factor for integration
  gp_names=unique(group)
  if(length(gp_names)!=2){stop("The number of groups must be 2.")}
  idx1 = group==gp_names[1]; idx2 = !idx1
  n1 = sum(idx1);n2=sum(idx2)
  
  if(any(!(rho_vec>0&rho_vec<=1))){stop("The threshold rho must be within zero and one.")}
  rho_names = paste0("rho=", rho_vec)
  names(rho_vec) = rho_names
  
  
  # ===============est_pre======
  # prelimiary steps: estimate mean/cov
  # centering by averages of each group
  X1bar = colMeans(X[idx1,]); X2bar = colMeans(X[idx2,])
  Xcent=X
  Xcent[idx1,]=t(t(X[idx1,])-X1bar)
  Xcent[idx2,]=t(t(X[idx2,])-X2bar)
  Y1bar = mean(Y[idx1]); Y2bar = mean(Y[idx2])
  Ybar_vec = Y; Ybar_vec[idx1] = Y1bar; Ybar_vec[idx2] = Y2bar
  Ycent = Y-Ybar_vec
  
  # pooled covariance estimation
  GaHat_pool = (cov(X1)*(n1-1)+cov(X2)*(n2-1)) / n
  eig_pool = eigen(GaHat_pool)
  gaHat_pool = eig_pool$values[1:h_max] * scal
  phiHat_pool = eig_pool$vectors[,1:h_max] / sqrt(scal)
  
  XprojCent_pool = Xcent %*% phiHat_pool * scal
  
  # separate estimations
  
  GaHat1 = cov(X1)*(n1-1)/n1; 
  # eig1 = eigen(GaHat1)
  GaHat2 = cov(X2)*(n2-1)/n2; 
  # eig2 = eigen(GaHat2)
  
  # for SimulDiag
  GaHat_sum = GaHat1 + GaHat2
  eig_sum = eigen(GaHat_sum)
  gaHat_sum = eig_sum$values[1:h_max] * scal
  phiHat_sum = eig_sum$vectors[,1:h_max] / sqrt(scal)
  
  GaHatInv_sum_sqrt_h = matrix(0, tt, tt)
  etaHatSD_list = list()
  gaHat_SD_list = phiHat_SD_list = XprojCentSD_list = list()
  for(h in 1:h_max){
    GaHatInv_sum_sqrt_each_j = outer(phiHat_sum[,h], phiHat_sum[,h])/sqrt(gaHat_sum[h])
    GaHatInv_sum_sqrt_h = GaHatInv_sum_sqrt_each_j + GaHatInv_sum_sqrt_h
    
    GaHatSD = GaHatInv_sum_sqrt_h %*% GaHat1 %*% GaHatInv_sum_sqrt_h * scal * scal
    
    # need to enforce symmetry when both ev and ef are unequal?
    GaHatSD[lower.tri(GaHatSD)] = t(GaHatSD)[lower.tri(GaHatSD)]
    
    # because ga_1 ga_2 are orthogonal, there are many ones and others are zero
    eigSD = eigen(GaHatSD)
    etaHatSD = eigSD$values[1:h] * scal
    # etaHatSD_list[[h]] = etaHatSD
    zetaHatSD = eigSD$vectors[,1:h] / sqrt(scal)
    
    gaHatSD = etaHatSD/(1-etaHatSD)
    phiHatSD = t(t(GaHatInv_sum_sqrt_h %*% zetaHatSD * scal) / sqrt(1-etaHatSD))
    
    XprojCentSD_h = Xcent %*% phiHatSD[,1:h] * scal
    
    # store
    etaHatSD_list[[h]] = etaHatSD
    gaHat_SD_list[[h]] = gaHatSD
    phiHat_SD_list[[h]] = phiHatSD
    XprojCentSD_list[[h]] = XprojCentSD_h
    
  }
  # need to use faster decaying eigenvalues for ga1
  # otherwise, the order of eta will be revsersed.
  
  
  # ============stat=======
  
  # test statistic for both equal/unequal X covariance but equal \e variance
  
  # for equal cov
  res_lm_pool_list =
    betaHatCoef_pool_list = list()
  MSE_pool_vec =
    TstatEq_vec = rep(0,h_max)
  Yhat_mat_pool = erHatCent_mat_pool = matrix(0, h_max, n)
  
  # for unequal cov
  res_lmSD_list =
    betaHatCoefSD_list = betaHatSD_list = list()
  MSE_SD_pool_vec =
    TstatUneq_vec = rep(0,h_max)
  Yhat_mat_SD = erHatCent_mat_SD = matrix(0, h_max, n)
  
  # for unequal cov, enforcing the null
  res_lmSD_pool_list =
    betaHatCoefSD_pool_list = list()
  betaHat_SD_pool_mat = matrix(0, h_max, tt)
  Yhat_SD_pool_mat = erHatCent_SD_pool_mat = matrix(0, h_max, n)
  
  # # extra statistics
  # TstatEq_Dt_vec = TstatEq_beta_vec = TstatEq_extreme_vec = rep(0,h_max)
  # TstatUneq_Dt_vec = TstatUneq_beta_vec = TstatUneq_extr_vec = rep(0,h_max)
  
  for(h in 1:h_max){
    
    # ========eq=========
    # test statistic for equal X covariance
    res_lm1_pool = lm(Ycent[idx1]~XprojCent_pool[idx1,1:h])
    res_lm2_pool = lm(Ycent[idx2]~XprojCent_pool[idx2,1:h])
    betaHatCeof1_pool_j = unname(coef(res_lm1_pool)[-1])
    betaHatCeof2_pool_j = unname(coef(res_lm2_pool)[-1])
    betaHatCoef_pool_diff = betaHatCeof1_pool_j-betaHatCeof2_pool_j
    
    SSE_pool = 
      summary(res_lm1_pool)$sigma^2*(n1-res_lm1_pool$rank)+
      summary(res_lm2_pool)$sigma^2*(n2-res_lm2_pool$rank)
    MSE_pool = SSE_pool / (n-h-1)
    SEsq_pool = (MSE_pool / gaHat_pool[1:h]) * (1/n1 + 1/n2)
    Tj_std_eq = 
      betaHatCoef_pool_diff / sqrt(SEsq_pool)
    TstatEq = sum(Tj_std_eq^2)
    
    # store
    # res_lm_pool_list[[h]] = list(gp1=res_lm1_pool, gp2=res_lm2_pool)
    betaHatCoef_pool_list[[h]] = rbind(
      gp1=betaHatCeof1_pool_j,
      gp2=betaHatCeof2_pool_j,
      diff = betaHatCoef_pool_diff,
      Tj_std_eq = Tj_std_eq
    )
    MSE_pool_vec[h] = MSE_pool
    Yhat_mat_pool[h,idx1] = unname(fitted(res_lm1_pool))
    Yhat_mat_pool[h,idx2] = unname(fitted(res_lm2_pool))
    erHatCent_mat_pool[h,idx1] = unname(resid(res_lm1_pool))
    erHatCent_mat_pool[h,idx2] = unname(resid(res_lm2_pool))
    TstatEq_vec[h] = TstatEq
    
    # # extra statistics
    # TstatEq_Dt_vec[h] = sum(gaHat_pool[1:h]*Tj_std_eq^2)
    # temp = Tj_std_eq^2/gaHat_pool[1:h]
    # TstatEq_beta_vec[h] = sum(temp)
    # TstatEq_extreme_vec[h] = sum(temp/gaHat_pool[1:h])
    
    
    # ==========uneq SD==========
    # uneq cov, using simultaneous diagonalization
    
    
    
    # 
    XprojCentSD_h = XprojCentSD_list[[h]]
    gaHatSD_h = gaHat_SD_list[[h]]
    phiHatSD_h = phiHat_SD_list[[h]]
    
    res_lm1SD = lm(Ycent[idx1]~XprojCentSD_h[idx1,1:h])
    res_lm2SD = lm(Ycent[idx2]~XprojCentSD_h[idx2,1:h])
    betaHatCeof1SD_j = unname(coef(res_lm1SD)[-1])
    betaHatCeof2SD_j = unname(coef(res_lm2SD)[-1])
    betaHatCoefSD_diff = betaHatCeof1SD_j-betaHatCeof2SD_j
    
    SSE_SD_pool =
      summary(res_lm1SD)$sigma^2*(n1-res_lm1SD$rank)+
      summary(res_lm2SD)$sigma^2*(n2-res_lm2SD$rank)
    MSE_SD_pool = SSE_SD_pool / (n-h-1)
    wtSD = ((n1*gaHatSD_h[1:h])^(-1) + (n2)^(-1))
    SEsq_SD_pool = MSE_SD_pool * wtSD
    Tj_std_uneq = betaHatCoefSD_diff / sqrt(SEsq_SD_pool)
    TstatUneq = sum(Tj_std_uneq^2)
    
    # store
    betaHatCoefSD_list[[h]] = rbind(
      gp1=betaHatCeof1SD_j,
      gp2=betaHatCeof2SD_j,
      diff = betaHatCoefSD_diff,
      Tj_std_uneq = Tj_std_uneq
    )
    MSE_SD_pool_vec[h] = MSE_SD_pool
    Yhat_mat_SD[h,idx1] = unname(fitted(res_lm1SD))
    Yhat_mat_SD[h,idx2] = unname(fitted(res_lm2SD))
    erHatCent_mat_SD[h,idx1] = unname(resid(res_lm1SD))
    erHatCent_mat_SD[h,idx2] = unname(resid(res_lm2SD))
    TstatUneq_vec[h] = TstatUneq
    
    # # extra statistics
    # TstatUneq_Dt_vec[h] = sum(wtSD*Tj_std_eq^2)
    # temp = Tj_std_eq^2/wtSD
    # TstatUneq_beta_vec[h] = sum(temp)
    # TstatUneq_extr_vec[h] = sum(temp/wtSD)
    
    
    # ====================
    
    # pooled estimates based on SD for enforcing the null
    res_lm_SD_pool = lm(Ycent~XprojCentSD_h)
    betaHatCeof_SD_pool_vec = unname(coef(res_lm_SD_pool)[-1])
    betaHat_SD_pool_h = 
      colSums(betaHatCeof_SD_pool_vec*t(phiHatSD_h[,1:h]))
    Yhat_SD_pool_h = 
      Ybar_vec + colSums(betaHatCeof_SD_pool_vec*t(XprojCentSD_h))
    erHatCent_SD_pool_h = Y-Yhat_SD_pool_h
    
    res_lmSD_pool_list[[h]] = res_lm_SD_pool
    betaHatCoefSD_pool_list[[h]] = betaHatCeof_SD_pool_vec
    betaHat_SD_pool_mat[h,] = betaHat_SD_pool_h
    Yhat_SD_pool_mat[h,] = Yhat_SD_pool_h
    erHatCent_SD_pool_mat[h,] = erHatCent_SD_pool_h
    
  }
  
  
  
  
  # ===============
  
  # choice of h_n using FVE like criteria
  
  
  gaHatPDsum_total = sum(eig_pool$values*scal)
  
  gaHatPDsum = cumsum(gaHat_pool)
  FVE_PD = (gaHatPDsum/gaHatPDsum_total)[1:h_max]
  RVE_PD = c(gaHatPDsum[1:(h_max-1)] / gaHatPDsum[2:h_max],0)
  
  gaHatSDsum = sapply(gaHat_SD_list, sum)
  RVE_SD = c(gaHatSDsum[1:(h_max-1)]/gaHatSDsum[2:h_max],0)
  
  
  # 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  Tstatftn_eq = function(Ycent_inner, h, diff_cent=0){
    res_lm1_pool = lm(Ycent_inner[idx1]~XprojCent_pool[idx1,1:h])
    res_lm2_pool = lm(Ycent_inner[idx2]~XprojCent_pool[idx2,1:h])
    betaHatCeof1_pool_j = unname(coef(res_lm1_pool)[-1])
    betaHatCeof2_pool_j = unname(coef(res_lm2_pool)[-1])
    betaHatCoef_pool_diff = betaHatCeof1_pool_j-betaHatCeof2_pool_j
    
    SSE_pool = 
      summary(res_lm1_pool)$sigma^2*(n1-res_lm1_pool$rank)+
      summary(res_lm2_pool)$sigma^2*(n2-res_lm2_pool$rank)
    MSE_pool = SSE_pool / (n-h-1)
    SEsq_pool = (MSE_pool / gaHat_pool[1:h]) * (1/n1 + 1/n2)
    Tj_std_eq = 
      (betaHatCoef_pool_diff-diff_cent) / sqrt(SEsq_pool)
    TstatEq = sum(Tj_std_eq^2)
    
    return(TstatEq)
  }
  Tstatftn_uneq = function(Ycent_inner, h, diff_cent=0){
    
    XprojCentSD_h = XprojCentSD_list[[h]]
    gaHatSD_h = gaHat_SD_list[[h]]
    
    res_lm1SD = lm(Ycent_inner[idx1]~XprojCentSD_h[idx1,1:h])
    res_lm2SD = lm(Ycent_inner[idx2]~XprojCentSD_h[idx2,1:h])
    betaHatCeof1SD_j = unname(coef(res_lm1SD)[-1])
    betaHatCeof2SD_j = unname(coef(res_lm2SD)[-1])
    betaHatCoefSD_diff = betaHatCeof1SD_j-betaHatCeof2SD_j
    
    SSE_SD_pool =
      summary(res_lm1SD)$sigma^2*(n1-res_lm1SD$rank)+
      summary(res_lm2SD)$sigma^2*(n2-res_lm2SD$rank)
    MSE_SD_pool = SSE_SD_pool / (n-h-1)
    wtSD = ((n1*gaHatSD_h[1:h])^(-1) + (n2)^(-1))
    SEsq_SD_pool = MSE_SD_pool * wtSD
    Tj_std_uneq = (betaHatCoefSD_diff-diff_cent) / sqrt(SEsq_SD_pool)
    TstatUneq = sum(Tj_std_uneq^2)
    
    return(TstatUneq)
  }
  
  
  
  
  # ===============================RB======
  # pooled regression estimator for residual bootstrap
  # # this might be appropriate only for equal cov?
  # for enforcing the null
  # using just pool looks better
  
  res_lm = lm(Ycent~XprojCent_pool)
  betaHatCeof_vec = unname(coef(res_lm)[-1])
  betaHat_mat = apply(betaHatCeof_vec*t(phiHat_pool), 2, cumsum)
  Yhat_pool_pool_mat = t(Ybar_vec + t(apply(betaHatCeof_vec*t(XprojCent_pool), 2, cumsum)))
  erHatCent_pool_pool_mat = t(Y-t(Yhat_pool_pool_mat))
  
  
  # residual bootstrap
  TstatStarEq_mat = TstatStarUneq_mat = 
    TstatStarEqH0_mat = TstatStarUneqH0_mat = 
    matrix(0, h_max, B)
  
  for(iBTS in 1:B){
    
    # ==========resample=========
    # resample bootstrap errors and responses
    
    idxBTS = ceiling(n*runif(n))
    
    # not enforcing the null: equal cov
    erStar_mat_eq = erHatCent_mat_pool[,idxBTS]
    
    YStar_mat_eq = erStar_mat_eq
    YStar_mat_eq[,idx1] = Yhat_mat_pool[,idx1] + erStar_mat_eq[,idx1]
    YStar_mat_eq[,idx2] = Yhat_mat_pool[,idx2] + erStar_mat_eq[,idx2]
    
    YStarCent_mat_eq = YStar_mat_eq
    YStarCent_mat_eq[,idx1] = YStar_mat_eq[,idx1] - rowMeans(YStar_mat_eq[,idx1])
    YStarCent_mat_eq[,idx2] = YStar_mat_eq[,idx2] - rowMeans(YStar_mat_eq[,idx2])
    
    # not enforcing the null: unequal cov
    erStar_mat_uneq = erHatCent_mat_SD[,idxBTS]
    
    YStar_mat_uneq = erStar_mat_uneq
    YStar_mat_uneq[,idx1] = Yhat_mat_SD[,idx1] + erStar_mat_uneq[,idx1]
    YStar_mat_uneq[,idx2] = Yhat_mat_SD[,idx2] + erStar_mat_uneq[,idx2]
    
    YStarCent_mat_uneq = YStar_mat_uneq
    YStarCent_mat_uneq[,idx1] = YStar_mat_uneq[,idx1] - rowMeans(YStar_mat_uneq[,idx1])
    YStarCent_mat_uneq[,idx2] = YStar_mat_uneq[,idx2] - rowMeans(YStar_mat_uneq[,idx2])
    
    # enforcing the null: equal cov
    # this might be appropriate only for equal cov? Seems not
    # enforceH0 seems to work for both
    
    erStar_mat_eq_H0= erHatCent_pool_pool_mat[,idxBTS]
    
    YStar_mat_eq_H0 = erStar_mat_eq_H0
    YStar_mat_eq_H0[,idx1] = Yhat_pool_pool_mat[,idx1] + erStar_mat_eq_H0[,idx1]
    YStar_mat_eq_H0[,idx2] = Yhat_pool_pool_mat[,idx2] + erStar_mat_eq_H0[,idx2]
    
    YStarCent_mat_eq_H0 = YStar_mat_eq_H0
    YStarCent_mat_eq_H0[,idx1] = YStar_mat_eq_H0[,idx1] - rowMeans(YStar_mat_eq_H0[,idx1])
    YStarCent_mat_eq_H0[,idx2] = YStar_mat_eq_H0[,idx2] - rowMeans(YStar_mat_eq_H0[,idx2])
    
    # enforcing the null: unequal cov
    
    erStar_mat_uneq_H0= erHatCent_SD_pool_mat[,idxBTS]
    
    YStar_mat_uneq_H0 = erStar_mat_uneq_H0
    YStar_mat_uneq_H0[,idx1] = Yhat_SD_pool_mat[,idx1] + erStar_mat_uneq_H0[,idx1]
    YStar_mat_uneq_H0[,idx2] = Yhat_SD_pool_mat[,idx2] + erStar_mat_uneq_H0[,idx2]
    
    YStarCent_mat_uneq_H0 = YStar_mat_uneq_H0
    YStarCent_mat_uneq_H0[,idx1] = YStar_mat_uneq_H0[,idx1] - rowMeans(YStar_mat_uneq_H0[,idx1])
    YStarCent_mat_uneq_H0[,idx2] = YStar_mat_uneq_H0[,idx2] - rowMeans(YStar_mat_uneq_H0[,idx2])
    
    # ======================
    
    for(h in 1:h_max){
      
      # not enforce H0, equal cov
      TstatStarEq_mat[h,iBTS] = 
        Tstatftn_eq(YStarCent_mat_eq[h,], h, betaHatCoef_pool_list[[h]]["diff",])
      # not enforce H0, unequal cov
      TstatStarUneq_mat[h,iBTS] = 
        Tstatftn_uneq(YStarCent_mat_uneq[h,], h, betaHatCoefSD_list[[h]]["diff",])
      # enforce H0, equal cov
      TstatStarEqH0_mat[h,iBTS] = 
        Tstatftn_eq(YStarCent_mat_eq_H0[h,], h)
      # enforce H0, unequal cov
      TstatStarUneqH0_mat[h,iBTS] = 
        Tstatftn_uneq(YStarCent_mat_uneq_H0[h,], h)
      
    }
    
  }
  
  # ====================================
  
  # ==================results===========
  stat = rbind(
    PD = TstatEq_vec  ,
    SD = TstatUneq_vec
  )
  colnames(stat) = paste0("h=",1:h_max)
  pval = abind::abind(
    PD = rbind(
      asymp = 1-pchisq(TstatEq_vec, 1:h_max),
      BTS = rowMeans(TstatEq_vec<TstatStarEq_mat),
      BTSH0 = rowMeans(TstatEq_vec<TstatStarEqH0_mat)
    ),
    SD = rbind(
      asymp = 1-pchisq(TstatUneq_vec, 1:h_max),
      BTS = rowMeans(TstatUneq_vec<TstatStarUneq_mat),
      BTSH0 = rowMeans(TstatUneq_vec<TstatStarUneqH0_mat)
    ),
    along=3
  )
  dimnames(pval)[[2]] = paste0("h=",1:h_max)
  VE = rbind(
    FVE_PD = FVE_PD,
    RVE_PD = RVE_PD,
    RVE_SD = RVE_SD
  )
  colnames(VE) = paste0("h=",1:h_max)
  
  
  
  hVEsel = function(VE, rho_vec){
    # VE: vector
    # rho_vec: vector
    # return: vector with length equal to length(rho_vec)
    sapply(rho_vec, function(rho){
      where = VE>rho;
      ifelse(any(where), min(which(where)), which.max(VE))
    })
  }
  hSel = t(apply(VE, 1, hVEsel, rho_vec=rho_vec))
  
  stat_hSel = rbind(
    FVE_PD = stat["PD",hSel["FVE_PD",]],
    RVE_PD = stat["PD",hSel["RVE_PD",]],
    RVE_SD = stat["SD",hSel["RVE_SD",]]
  )
  colnames(stat_hSel) = rho_names
  pval_hSel = abind::abind(
    FVE_PD = pval[,hSel["FVE_PD",],"PD"],
    RVE_PD = pval[,hSel["RVE_PD",],"PD"],
    RVE_SD = pval[,hSel["RVE_SD",],"SD"],
    along=3
  )
  dimnames(pval_hSel)[[2]] = rho_names
  
  
  res_hEach = list( stat = stat, pval=pval, VE = VE )
  res_hSel = list( stat = stat_hSel, pval = pval_hSel, hSel = hSel)
  # ===================================
  
  res = list(res_hEach = res_hEach, res_hSel = res_hSel)
  res
  
}