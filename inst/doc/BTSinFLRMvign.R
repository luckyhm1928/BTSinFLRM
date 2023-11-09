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

## ---- echo=FALSE, out.width='60%', fig.align='center', fig.width=7, fig.height=7----
matplot(tGrid, t(X), type='l', xlab="", ylab="",
        main="Generated Gaussian regressor curves")

## -----------------------------------------------------------------------------
# trunction parameters
kn_vec = c(3); hn_vec = c(4,5,6); gn_vec = c(4,5)

# residual bootstrap
M_bts=100

library(BTSinFLRM)

## ---- echo=FALSE, out.width='60%', fig.align='center', fig.width=7, fig.height=7----
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

