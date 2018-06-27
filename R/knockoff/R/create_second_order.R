#' Second-order Gaussian knockoffs
#' 
#' This function samples second-order multivariate Gaussian knockoff variables.
#' First, a multivariate Gaussian distribution is fitted to the observations of X.
#' Then, Gaussian knockoffs are generated according to the estimated model.
#' 
#' @param X n-by-p matrix of original variables.
#' @param method either "equi", "sdp" or "asdp" (default: "asdp").
#' This determines the method that will be used to minimize the correlation between the original variables and the knockoffs.
#' @param shrink whether to shrink the estimated covariance matrix (default: F).
#' @return A n-by-p matrix of knockoff variables.
#'  
#' @family create
#' 
#' @details
#' If the argument \code{shrink} is set to T, a James-Stein-type shrinkage estimator for
#' the covariance matrix is used instead of the traditional maximum-likelihood estimate. This option
#' requires the package \code{corpcor}. See \code{\link[corpcor]{cov.shrink}} for more details.
#' 
#' Even if the argument \code{shrink} is set to F, in the case that the estimated covariance 
#' matrix is not positive-definite, this function will apply some shrinkage.
#' 
#' @references 
#'   Candes et al., Panning for Gold: Model-free Knockoffs for High-dimensional Controlled Variable Selection,
#'   arXiv:1610.02351 (2016).
#'   \href{https://web.stanford.edu/group/candes/knockoffs/index.html}{https://web.stanford.edu/group/candes/knockoffs/index.html}
#'   
#' @examples
#' p=200; n=100; k=15
#' rho = 0.4
#' Sigma = toeplitz(rho^(0:(p-1)))
#' X = matrix(rnorm(n*p),n) %*% chol(Sigma)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#' 
#' # Basic usage with default arguments
#' result = knockoff.filter(X, y, knockoffs=create.second_order)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' knockoffs = function(X) create.second_order(X, method='equi')
#' result = knockoff.filter(X, y, knockoffs=knockoffs)
#' print(result$selected)   
#'   
#' @export
create.second_order <- function(X, method=c("asdp","equi","sdp"), shrink=F) {
  method = match.arg(method)
  
  # Estimate the mean vectorand covariance matrix
  mu = colMeans(X)
  
  # Estimate the covariance matrix
  if(!shrink) {
    Sigma = cov(X)
    # Verify that the covariance matrix is positive-definite
    if(!is_posdef(Sigma)) {
      shrink=TRUE
    }
  }
  if(shrink) {
    if (!requireNamespace('corpcor', quietly=T))
      stop('corpcor is not installed', call.=F)
    Sigma = tryCatch({suppressWarnings(matrix(as.numeric(corpcor::cov.shrink(X,verbose=F)), nrow=ncol(X)))},
                     warning = function(w){}, error = function(e) {
                       stop("SVD failed in the shrinkage estimation of the covariance matrix. Try upgrading R to version >= 3.3.0")
                     }, finally = {})
  }

  # Sample the Gaussian knockoffs
  create.gaussian(X, mu, Sigma, method=method)
}