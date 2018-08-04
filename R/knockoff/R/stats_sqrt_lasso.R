#' Importance statistics based on the square-root lasso
#' 
#' Computes the signed maximum statistic
#'   \deqn{W_j = \max(Z_j, \tilde{Z}_j) \cdot \mathrm{sgn}(Z_j - \tilde{Z}_j),}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the maximum values of 
#' \eqn{\lambda} at which the jth variable and its knockoff, respectively,
#' enter the SQRT lasso model.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables of numeric type.
#' @param ... additional arguments specific to \code{slim}.
#' @return A vector of statistics \eqn{W} of length p.
#' 
#' @details With default parameters, this function uses the package \code{flare}
#' to run the SQRT lasso. By specifying the appropriate optional parameters, 
#' one can use different Lasso variants including Dantzig Selector, LAD Lasso,
#' SQRT Lasso and Lq Lasso for estimating high dimensional sparse linear models.
#' 
#' For a complete list of the available additional arguments, see \code{\link[flare]{slim}}.
#' 
#' @family statistics
#' 
#' @examples
#' p=50; n=50; k=10
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#' knockoffs = function(X) create.gaussian(X, mu, Sigma)
#' 
#' # Basic usage with default arguments
#' result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=stat.sqrt_lasso)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = stat.sqrt_lasso
#' k_stat = function(X, X_k, y) foo(X, X_k, y, q=0.5)
#' result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname stat.sqrt_lasso
#' @export
stat.sqrt_lasso <- function(X, X_k, y, ...) {
  if (!requireNamespace('flare', quietly=T))
    stop('flare is not installed', call.=F)
  if (!(is.vector(y) && is.numeric(y)))  {
    stop('Knockoff statistic stat.sqrt_lasso requires the input y to be a numeric vector')
  }
  p = ncol(X)
  
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)
  
  # Compute statistics
  Z = sqrt_lasso_coeffs(cbind(X.swap, Xk.swap), y, ...)
  p = ncol(X)
  orig = 1:p
  W = pmax(Z[orig], Z[orig+p]) * sign(Z[orig] - Z[orig+p])
  
  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
}

#' SQRT lasso
#' 
#' Perform variable selection with SQRT lasso
#' 
#' @param X matrix of predictors
#' @param y response vector
#' @return vector with jth component the selection probability of variable j
#' 
#' @keywords internal
sqrt_lasso_coeffs <- function(X, y, nlambda=NULL, ...) {
  X = scale(X)
  n = nrow(X)
  fit = flare::slim(X, y)
  first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
  indices <- apply(fit$beta, 1, first_nonzero)
  names(indices) <- NULL
  ifelse(is.na(indices), 0, fit$lambda[indices] * n)
}
