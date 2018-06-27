#' Importance statistics based the lasso with cross-validation
#' 
#' Fits a linear regression model via penalized maximum likelihood and cross-validation.
#' Then, compute the difference statistic
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the coefficient estimates for the 
#' jth variable and its knockoff, respectively. The value of the regularization
#' parameter \eqn{\lambda} is selected by cross-validation and computed with \code{glmnet}.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables. It should be numeric
#' @param cores Number of cores used to compute the statistics by running cv.glmnet.
#' If not specified, the number of cores is set to approximately half of the number of cores 
#' detected by the parallel package.
#' @param ... additional arguments specific to \code{glmnet} (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#' 
#' @details This function uses the \code{glmnet} package to fit the lasso path and 
#' is a wrapper around the more general \link{stat.glmnet_coefdiff}.
#' 
#' The statistics \eqn{W_j} are constructed by taking the difference 
#' between the coefficient of the j-th variable and its knockoff.
#'  
#' By default, the value of the regularization parameter is chosen by 10-fold cross-validation.
#' 
#' The optional \code{nlambda} parameter can be used to control the granularity of the 
#' grid of \eqn{\lambda}'s. The default value of \code{nlambda} is \code{500},
#' where \code{p} is the number of columns of \code{X}.
#' 
#' Unless a lambda sequence is provided by the user, this function generates it on a 
#' log-linear scale before calling 'glmnet' (default 'nlambda': 500).
#' 
#' For a complete list of the available additional arguments, see \code{\link[glmnet]{cv.glmnet}}
#' and \code{\link[glmnet]{glmnet}}.
#' 
#' @family statistics
#' 
#' @examples
#' p=200; n=100; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#' knockoffs = function(X) create.gaussian(X, mu, Sigma)
#' 
#' # Basic usage with default arguments
#' result = knockoff.filter(X, y, knockoffs=knockoffs, 
#'                            statistic=stat.lasso_coefdiff)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = stat.lasso_coefdiff
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname stat.lasso_coefdiff
#' @export
stat.lasso_coefdiff <- function(X, X_k, y, cores=2, ...) {
  if( is.numeric(y) ){
    y = as.vector(y)
  } else {
    stop('Knockoff statistic stat.lasso_coefdiff requires the input y to be a numeric vector')
  }
  
  stat.glmnet_coefdiff(X, X_k, y, family='gaussian', cores=cores, ...)
}
