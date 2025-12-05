#' @title Bootstrap inference in functional linear regression models with scalar response.
#' 
#' @author Hyemin Yeon \email{hyeon@@kent.edu}
#' 
#' @description
#' 
#' The package implements the residual, paired, and (mutliplier) wild bootstrap methods described in Yeon, Dai, and Nordman (2023, 2024, 2025);
#' these are developed for inference on projections of the slope function onto new predictors/regressors 
#' (i.e., the mean responses at the new regressors) 
#' in functional linear regression models (FLRMs) with scalar response.
#' The residual bootstrap is designed only for homoscedastic FLRMs
#' while the paired and wild bootstraps are developed to deal with heteroscedastic error cases.
#' \code{RBinFLRM} provides confidence and prediction intervals for projections 
#' based on residual bootstrap proposed by Yeon, Dai, and Nordman (2023).
#' \code{PBinFLRM} performs a modified paired bootstrap method by Yeon, Dai, and Nordman (2024)
#' to find confidence intervals for projections
#' and to compute p-values from the hypothesis tests 
#' of whether the projections onto multiple new regressors are simultaneously zero or not.
#' \code{WBinFLRM} shows the same results as the ones from \code{PBinFLRM} 
#' but based on (multiplier) wild bootstrap developed in Yeon, Dai, and Nordman (2025).
#' For more details see the help vignette:
#' \code{vignette("BTSinFLRMvign", package = "BTSinFLRM")}.
#' 
#' 
#' @docType package
#' @name BTSinFLRM
#' @useDynLib BTSinFLRM
#' @import Rcpp
#' 
#' 
#' 
#' @seealso
#' \code{\link{RBinFLRM}}
#' \code{\link{PBinFLRM}}
#' \code{\link{WBinFLRM}}
#' 
#' @examples 
#' vignette("BTSinFLRMvign", package = "BTSinFLRM")
#' 
#' @references 
#' 
#' Yeon, H., Dai, X., and Nordman, D. (2023). Bootstrap inference in functional linear regression models. Bernoulli
#' 
#' Yeon, H., Dai, X., and Nordman, D. (2024). Bootstrap inference in functional linear regression models under heteroscedasticity. Electronic Journal of Statistics
#' 
#' Yeon, H., Dai, X., and Nordman, D. (2025). Different bootstrap methods in functional linear regression models. Submitted
#' 
NULL
