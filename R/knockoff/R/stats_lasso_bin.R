#' Importance statistics based on regularized logistic regression
#' 
#' Fit the lasso path and computes the difference statistic
#'   \deqn{W_j = Z_j - \tilde{Z}_j}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of the 
#' regularization parameter \eqn{\lambda} at which the jth variable 
#' and its knockoff enter the penalized logistic regression model, respectively.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables. It should be either a factor with two levels, 
#' or a two-column matrix of counts or proportions 
#' (the second column is treated as the target class; for a factor, the last level 
#' in alphabetical order is the target class). If y is presented as a vector, 
#' it will be coerced into a factor.
#' @param ... additional arguments specific to \code{glmnet} (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#' 
#' @details This function uses \code{glmnet} to compute the lasso path
#' on a fine grid of \eqn{\lambda}'s.
#' 
#' The \code{nlambda} parameter can be used to control the granularity of the 
#' grid of \eqn{\lambda}'s. The default value of \code{nlambda} is \code{500}.
#' 
#' This function is a wrapper around the more general \code{\link{stat.glmnet_lambdadiff}}.
#' 
#' For a complete list of the available additional arguments, see \code{\link[glmnet]{glmnet}}
#' or \code{\link[lars]{lars}}.
#' 
#' @family statistics
#' 
#' @examples
#' set.seed(2022)
#' p=200; n=100; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' pr = 1/(1+exp(-X %*% beta))
#' y = rbinom(n,1,pr)
#' knockoffs = function(X) create.gaussian(X, mu, Sigma)
#' 
#' # Basic usage with default arguments
#' result = knockoff.filter(X, y, knockoffs=knockoffs, 
#'                            statistic=stat.lasso_lambdadiff_bin)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = stat.lasso_lambdadiff_bin
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname stat.lasso_lambdadiff_bin
#' @export
stat.lasso_lambdadiff_bin <- function(X, X_k, y, ...) {
  stat.glmnet_lambdadiff(X, X_k, y, family='binomial', ...)
}

#' Penalized logistic regression statistics for knockoff
#' 
#' Computes the signed maximum statistic
#'   \deqn{W_j = \max(Z_j, \tilde{Z}_j) \cdot \mathrm{sgn}(Z_j - \tilde{Z}_j),}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of 
#' \eqn{\lambda} at which the jth variable and its knockoff, respectively,
#' enter the penalized logistic regression model.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables. It should be either a factor with two levels, 
#' or a two-column matrix of counts or proportions 
#' (the second column is treated as the target class; for a factor, the last level 
#' in alphabetical order is the target class). If y is presented as a vector, 
#' it will be coerced into a factor.
#' @param ... additional arguments specific to \code{glmnet} or \code{lars} (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#'   
#' @details This function uses \code{glmnet} to compute the regularization path
#' on a fine grid of \eqn{\lambda}'s.
#' 
#' The additional \code{nlambda} 
#' parameter can be used to control the granularity of the grid of \eqn{\lambda} values. 
#' The default value of \code{nlambda} is \code{500}.
#' 
#' This function is a wrapper around the more general 
#' \link{stat.glmnet_lambdadiff}.
#' 
#' For a complete list of the available additional arguments, see \code{\link[glmnet]{glmnet}}.
#' 
#' @examples
#' p=200; n=100; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' pr = 1/(1+exp(-X %*% beta))
#' y = rbinom(n,1,pr)
#' knockoffs = function(X) create.gaussian(X, mu, Sigma)
#' 
#' # Basic usage with default arguments
#' result = knockoff.filter(X, y, knockoff=knockoffs,
#'                            statistic=stat.lasso_lambdasmax_bin)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = stat.lasso_lambdasmax_bin
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname stat.lasso_lambdasmax_bin
#' @export
stat.lasso_lambdasmax_bin <- function(X, X_k, y, ...) {
  stat.glmnet_lambdasmax(X, X_k, y, family='binomial', ...)
}