% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PBinFLRM.R
\name{PBinFLRM}
\alias{PBinFLRM}
\title{Paired bootstrap in functional linear regression models}
\usage{
PBinFLRM(X, Y, X0, tGrid, kn_vec, hn_vec, gn_vec, M_bts = 1000, ap = 0.05)
}
\arguments{
\item{X}{An n by p matrix of regressor curves. Each row represents one observed regressor curve.}

\item{Y}{An vector of n responses.}

\item{X0}{An n0 by p matrix of new regressor curves. Each row represents one new regressor curve.}

\item{tGrid}{A vector of p densely equi-spaced time grid points where the curve are evaluated.}

\item{kn_vec}{A vector of truncation level k associated with the estimated error variance.}

\item{hn_vec}{A vector of truncation level h associated with estimation of the slope function.}

\item{gn_vec}{A vector of truncation level g where the estimated slope function plays a role of the true parameter in the bootstrap world.}

\item{M_bts}{An integer. The number of bootstrap resamples with default 1000.}

\item{ap}{A numeric scalar. The significance level for inference.}
}
\value{
A list containing the following fields:
\item{est}{The FPCR estimated slope function with different trunction levels} 
\item{CLT}{The individual intervals based on CLT with different trunction levels}   
\item{PBintervals}{The intervals based on paired bootstrap with different truncation levels. 
Intervals are provided with or without either symmetrization or studentization.}
\item{PBpvalues}{P-values from the bootstrap hypothesis tests 
of whether the projections onto multiple new regressors are simultaneously zero or not.}
}
\description{
Perform paired bootstrap inference for projections of new regressors onto the slope function 
in functional linear regression models with scalar response under heteroscedasticity.

Apply the central limit theorem (CLT) (Yeon, Dai, and Nordman, 2024a) to find individual confidence intervals.
Apply the paired bootstrap (Yeon, Dai, and Nordman, 2024a) to find individual and simultaneous confidence intervals.
Intervals are provided with or without either symmetrization or studentization.
Simultaneous intervals are studentized by either the estimated scaling or the bootstrap scaling.
Apply the bootstrap hypothesis tests (Yeon, Dai, and Nordman, 2024a) 
of whether the projections onto multiple new regressors are simultaneously zero or not.
}
\examples{

# example

library(BTSinFLRM)

# simulate the regressor curves from Karhunen-Loeve expansion and responses
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
beta = c(bj \%*\% phi)

# generate (Gaussian) regressors and responses
n=50; n0=3
xi = matrix(rnorm(n*J), n, J); X = xi \%*\% (phi * sqrt(ga))
xi0 = matrix(rnorm(n0*J), n0, J); X0 = xi0 \%*\% (phi * sqrt(ga))
er = runif(n); Y = c(X \%*\% beta / tt + er)

# trunction parameters
kn_vec = c(3,4); hn_vec = c(4,5,6); gn_vec = c(5)

# paired bootstrap
M_bts=100
PBinFLRM(X, Y, X0, tGrid, kn_vec, hn_vec, gn_vec, M_bts)

}
\seealso{
\code{\link{RBinFLRM}}
\code{\link{WBinFLRM}}
}
