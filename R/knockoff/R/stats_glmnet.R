#' Importance statistics based on a GLM
#' 
#' Fits a generalized linear model via penalized maximum likelihood and
#' computes the difference statistic
#'   \deqn{W_j = Z_j - \tilde{Z}_j}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of the 
#' regularization parameter \eqn{\lambda} at which the jth variable 
#' and its knockoff enter the model, respectively.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables. Quantitative for family="gaussian", 
#' or family="poisson" (non-negative counts). For family="binomial" 
#' should be either a factor with two levels, or a two-column matrix of counts 
#' or proportions (the second column is treated as the target class; for a factor, 
#' the last level in alphabetical order is the target class). For family="multinomial", 
#' can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. 
#' For either "binomial" or "multinomial", if y is presented as a vector, it will 
#' be coerced into a factor. For family="cox", y should be a two-column matrix with 
#' columns named 'time' and 'status'. The latter is a binary variable, with '1' 
#' indicating death, and '0' indicating right censored. The function Surv() in 
#' package survival produces such a matrix. For family="mgaussian", y is a matrix 
#' of quantitative responses.
#' @param family response type (see above).
#' @param ... additional arguments specific to \code{glmnet} (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#' 
#' @details This function uses \code{glmnet} to compute the regularization path
#' on a fine grid of \eqn{\lambda}'s.
#' 
#' The \code{nlambda} parameter can be used to control the granularity of the 
#' grid of \eqn{\lambda}'s. The default value of \code{nlambda} is \code{500}.
#' 
#' If the family is 'binomial' and a lambda sequence is not provided by the user, 
#' this function generates it on a log-linear scale before calling 'glmnet'.
#' 
#' The default response family is 'gaussian', for a linear regression model.
#' Different response families (e.g. 'binomial') can be specified by passing an
#' optional parameter 'family'.
#' 
#' For a complete list of the available additional arguments, see \code{\link[glmnet]{glmnet}}.
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
#'                            statistic=stat.glmnet_lambdadiff)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = stat.glmnet_lambdadiff
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname stat.glmnet_lambdadiff
#' @export
stat.glmnet_lambdadiff <- function(X, X_k, y, family='gaussian', ...) {
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)
  
  # Compute statistics
  Z = lasso_max_lambda(cbind(X.swap, Xk.swap), y, method='glmnet', family=family, ...)
  p = ncol(X)
  orig = 1:p
  W = Z[orig] - Z[orig+p]
  
  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
}

#' GLM statistics for knockoff
#' 
#' Computes the signed maximum statistic
#'   \deqn{W_j = \max(Z_j, \tilde{Z}_j) \cdot \mathrm{sgn}(Z_j - \tilde{Z}_j),}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of 
#' \eqn{\lambda} at which the jth variable and its knockoff, respectively,
#' enter the generalized linear model.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables. Quantitative for family="gaussian", 
#' or family="poisson" (non-negative counts). For family="binomial" 
#' should be either a factor with two levels, or a two-column matrix of counts 
#' or proportions (the second column is treated as the target class; for a factor, 
#' the last level in alphabetical order is the target class). For family="multinomial", 
#' can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. 
#' For either "binomial" or "multinomial", if y is presented as a vector, it will 
#' be coerced into a factor. For family="cox", y should be a two-column matrix with 
#' columns named 'time' and 'status'. The latter is a binary variable, with '1' 
#' indicating death, and '0' indicating right censored. The function Surv() in 
#' package survival produces such a matrix. For family="mgaussian", y is a matrix 
#' of quantitative responses.
#' @param family response type (see above).
#' @param ... additional arguments specific to \code{glmnet} (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#'   
#' @details This function uses \code{glmnet} to compute the regularization path
#' on a fine grid of \eqn{\lambda}'s.
#' 
#' The additional \code{nlambda} 
#' parameter can be used to control the granularity of the grid of \eqn{\lambda} values. 
#' The default value of \code{nlambda} is \code{500}.
#' 
#' If the family is 'binomial' and a lambda sequence is not provided by the user, 
#' this function generates it on a log-linear scale before calling 'glmnet'.
#' 
#' For a complete list of the available additional arguments, see \code{\link[glmnet]{glmnet}}.
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
#' result = knockoff.filter(X, y, knockoff=knockoffs,
#'                            statistic=stat.glmnet_lambdasmax)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = stat.glmnet_lambdasmax
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname stat.glmnet_lambdasmax
#' @export
stat.glmnet_lambdasmax <- function(X, X_k, y, family='gaussian', ...) {
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)
  
  # Compute statistics
  Z = lasso_max_lambda(cbind(X.swap, Xk.swap), y, method='glmnet', family=family, ...)
  p = ncol(X)
  orig = 1:p
  W = pmax(Z[orig], Z[orig+p]) * sign(Z[orig] - Z[orig+p])
  
  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
}

#' @keywords internal
lasso_max_lambda_lars <- function(X, y, ...) {
  if (!requireNamespace('lars', quietly=T))
    stop('lars is not installed', call.=F)
  
  fit <- lars::lars(X, y, normalize=T, intercept=F, ...)
  lambda <- rep(0, ncol(X))
  for (j in 1:ncol(X)) {
    entry <- fit$entry[j]
    if (entry > 0) lambda[j] <- fit$lambda[entry]
  }
  return(lambda)
}

#' @keywords internal
lasso_max_lambda_glmnet <- function(X, y, nlambda=500, intercept=T, standardize=T, ...) {
  if (!requireNamespace('glmnet', quietly=T))
    stop('glmnet is not installed', call.=F)
  
  # Standardize the variables
  if( standardize ){
    X = scale(X)
  }
    
  n = nrow(X); p = ncol(X)
  if (!methods::hasArg(family) ) family = "gaussian"
  else family = list(...)$family
  
  if (!methods::hasArg(lambda) ) {
    if( identical(family, "gaussian") ) {
      if(!is.numeric(y)) {
        stop('Input y must be numeric.')
      }
      # Unless a lambda sequence is provided by the user, generate it
      lambda_max = max(abs(t(X) %*% y)) / n
      lambda_min = lambda_max / 2e3
      k = (0:(nlambda-1)) / nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }

  fit <- glmnet::glmnet(X, y, lambda=lambda, intercept=intercept, 
                        standardize=F, standardize.response=F, ...)
  
  first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
  indices <- apply(fit$beta, 1, first_nonzero)
  names(indices) <- NULL
  ifelse(is.na(indices), 0, fit$lambda[indices] * n)
}

#' Maximum lambda in lasso model
#' 
#' Computes the earliest (largest) lambda's for which predictors enter the
#' lasso model.
#' 
#' @param X matrix of predictors
#' @param y response vector
#' @param method either 'glmnet' or 'lars'
#' @return vector of maximum lambda's
#' 
#' @keywords internal
lasso_max_lambda <- function(X, y, method=c('glmnet','lars'), ...) {
  switch(match.arg(method), 
         glmnet = lasso_max_lambda_glmnet(X,y,...),
         lars = lasso_max_lambda_lars(X,y,...)
         )
}