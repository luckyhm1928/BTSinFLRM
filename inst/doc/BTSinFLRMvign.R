## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----genData------------------------------------------------------------------
# simulate the regressor curves from Karhunen–Loève expansion and responses
J=15

# eigenvalues
c=2; a=3.5; dt = c*(1:J)^(-a)
g1 = c*VGAM::zeta(a)
ga = c(g1, sapply(1:(J-1), function(j){g1 - sum(dt[1:j])}))

# eigenfunctions
tt = 100; tGrid = seq(0, 1, len=tt); J=15
phi= t(fda::fourier(tGrid, J))

# slope function
b=3.5; bj = c*((1:J)^(-b))*(rbinom(J, 1, 0.5)*2-1)
beta = c(bj %*% phi)

# generate (Gaussian) regressors and responses
n=50; n0=3
xi = matrix(rnorm(n*J), n, J); X = xi %*% (phi * sqrt(ga))
xi0 = matrix(rnorm(n0*J), n0, J); X0 = xi0 %*% (phi * sqrt(ga))
er = runif(n); Y = c(X %*% beta / tt + er)

## ----echo=FALSE, out.width='60%', fig.align='center', fig.width=7, fig.height=7----
matplot(tGrid, t(X), type='l', xlab="", ylab="",
        main="Generated Gaussian regressor curves")

## -----------------------------------------------------------------------------
# trunction parameters
kn_vec = c(3); hn_vec = c(4,5,6); gn_vec = c(4,5)

# residual bootstrap
M_bts=100

library(BTSinFLRM)

## ----echo=FALSE, out.width='60%', fig.align='center', fig.width=7, fig.height=7----
resInfer = BTSinFLRM:::inferFLRM(X, Y, X0, tGrid, kn_vec, hn_vec, gn_vec)
tmx = max(c(kn_vec, hn_vec, gn_vec))
matplot(tGrid, t(resInfer$est$betaHat), type='l', xlab="", ylab="", 
        main="Estimated slope functions")
legend("bottomright", legend = paste0("kn=",1:tmx),
       lty=1:tmx, col=1:tmx, lwd=1)

## -----------------------------------------------------------------------------
resRB = RBinFLRM(X, Y, X0, tGrid, kn_vec, hn_vec, gn_vec, M_bts)

## -----------------------------------------------------------------------------
resRB$CLT

## -----------------------------------------------------------------------------
resRB$RB
# Access a component of list:
resRB$RB[["kn=3", "hn=5", "gn=4"]]

## -----------------------------------------------------------------------------
resPB = PBinFLRM(X, Y, X0, tGrid, kn_vec, hn_vec, gn_vec, M_bts)

## -----------------------------------------------------------------------------
resPB$CLT

## -----------------------------------------------------------------------------
resPB$PBintervals
# Access a component of list:
resPB$PBintervals[["kn=3", "hn=5", "gn=4"]]

## -----------------------------------------------------------------------------
resPB$PBpvalues
# Access a component of list:
resPB$PBpvalues[["kn=3", "hn=5", "gn=4"]]

## -----------------------------------------------------------------------------
resWB = WBinFLRM(X, Y, X0, tGrid, kn_vec, hn_vec, gn_vec, M_bts, multiplier="normal")

## -----------------------------------------------------------------------------
resWB$WBintervals
# Access a component of list:
resWB$WBintervals[["kn=3", "hn=5", "gn=4"]]

## -----------------------------------------------------------------------------
resWB$WBpvalues
# Access a component of list:
resWB$WBpvalues[["kn=3", "hn=5", "gn=4"]]

## ----test2sofr_evef-----------------------------------------------------------
Jtrue=20 # The number of basis elements used in data generation

# create eigenvalues based on polynomial eigengaps with rate a
ev_poly = function(a, c, J){
  dt = c*(1:J)^(-a)
  ld1 = c*VGAM::zeta(a)
  ld = c(ld1, sapply(1:(J-1), function(j){ld1 - sum(dt[1:j])}))
  return(ld)
}

# eigenvalues for unequal cov
ga_rough  = ev_poly(2.5, 2, Jtrue) # a=2.5
ga_smooth = ev_poly(5, 2, Jtrue)   # a=5

# Equi-spaced time grid points in [0,1]
tt = 50; tGrid_ttt = seq(0, 1, len=tt+1);
tGrid = tGrid_ttt[1:tt]
diffrange = diff(tGrid)[1] + diff(range(tGrid)) # when using left points
scal = diffrange / tt   # scaling factor for integration

# basis functions for beta, and for X under equal cov
phi_base = t(fda::fourier(tGrid_ttt, 2*Jtrue))[-1,1:tt]  # trigonometric functions
phi = phi_base[1:Jtrue,]

# basis functions for X under unequal cov
mono = t(sapply(1:Jtrue, function(j)tGrid^j))
phi_mono = BTSinFLRM:::orthoL2Equidense(mono, tGrid)
cheb = t(sapply(1:Jtrue, function(j)pracma::chebPoly(j, 2*tGrid-1)))
phi_cheb = BTSinFLRM:::orthoL2Equidense(cheb, tGrid)

# Determine the second order structures
eveq = TRUE # eigenvalues are equal?
efeq = TRUE # eigenfunctions are equal?
if(eveq){  # equal eigenvalues
  ga1 = ga2 = ga_rough
}else{  # unequal eigenvalues
  ga1 = ga_smooth; ga2 = ga_rough
}
if(efeq){  # equal eigenfunctions
  phi1 = phi; phi2 = phi
}else{  # unequal eigenfunctions
  phi1 = phi_mono; phi2 = phi_cheb
}

## -----------------------------------------------------------------------------
# Slope functions
b=2; Wbetaj = rbinom(Jtrue, 1, 0.5)*2-1
bj = 3*((1:Jtrue)^(-b))*Wbetaj; beta = c(bj %*% phi)

beta_arr = array(0, dim=c(2, tt, 2), dimnames=list(
  1:2, NULL, paste0("type", 0:1)
))
# 0: linear difference
beta_arr[1,,"type0"] = 1*beta
beta_arr[2,,"type0"] = 2*beta

# 1: orthogonal difference
beta_arr[1,,"type1"] = c(bj[1:5] %*% phi[1:5,])
beta_arr[2,,"type1"] = c(bj[1:5] %*% phi[6:10,])

beta1 = beta_arr[1,,"type0"]
c_alter = 0.5 # strength of alternative
beta2 = (1-c_alter)*beta1 + c_alter*beta_arr[2,,"type0"]

# Other factors

n=50; n1 = n2 = n/2      # n1=48; n2=52
group = c(rep(1,n1), rep(0,n2))
idx1 = as.logical(group); idx2 = !idx1
xi_type=0 # Gaussian random function
er_type=0 # Normal error
sig_sq =1 # the error variance

# # # # # # # # # # #
# generate data

# FPC scores
Wj = matrix(rnorm(n*Jtrue), n, Jtrue)
if(xi_type==2){
  xi = rexp(n) - 1 # normal*exp FPC scores
}else if(xi_type==1){
  xi = rnorm(n) # normal*normal FPC scores
}else{
  xi = rep(1,n) # independent FPC scores
}

# # unequal FPC scores
# xi1 = rnorm(n1); xi2 = rexp(n2) - 1
# xi[idx1,] = xi1; xi[idx2,] = xi2

Xfpc = matrix(0, n, Jtrue)
Xfpc[idx1,] = t(sqrt(ga1) * t(xi[idx1]*Wj[idx1,]))
Xfpc[idx2,] = t(sqrt(ga2) * t(xi[idx2]*Wj[idx2,]))

# regerssor functions

X = matrix(0, n, tt)
X[idx1,] = Xfpc[idx1,] %*% phi1
X[idx2,] = Xfpc[idx2,] %*% phi2

X1 = X[idx1,]; X2 = X[idx2,]

X1beta1 = c(X1 %*% beta1 * scal)
X2beta2 = c(X2 %*% beta2 * scal)

if(er_type==1){   # centered exponential error
  par_scale = sqrt(sig_sq)
  er = rexp(n, 1/par_scale) - par_scale
}else{     # normal error
  er = rnorm(n, 0,sd=sqrt(sig_sq))
}

# response variables
Y1 = X1beta1 + er[idx1]
Y2 = X2beta2 + er[idx2]
Y = rep(0,n); Y[idx1] = Y1; Y[idx2] = Y2

# Bootstrap tests for two-sample scalar-on-function regression (not run)
res2sofr = test2sofr(Y, X, tGrid, group, B=100)

## ----test2sofr_results--------------------------------------------------------
res2sofr$res_hSel$pval

