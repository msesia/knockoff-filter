#' Importance statistics based on stability selection
#' 
#' Computes the difference statistic
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are measure the importance
#' of the jth variable and its knockoff, respectively, based on the 
#' stability of their selection upon subsampling of the data.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y response vector (length n)
#' @param fitfun fitfun a function that takes the arguments x, y as above, 
#' and additionally the number of variables to include in each model q. 
#' The function then needs to fit the model and to return a logical vector 
#' that indicates which variable was selected (among the q selected variables).
#' The name of the function should be prefixed by 'stabs::'.
#' @param ... additional arguments specific to 'stabs' (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#'   
#' @details This function uses the \code{stabs} package to compute
#' variable selection stability. The selection stability of the j-th 
#' variable is defined as its probability of being selected upon random
#' subsampling of the data. The default method for selecting variables 
#' in each subsampled dataset is \code{\link[stabs]{glmnet.lasso_maxCoef}}.
#' 
#' For a complete list of the available additional arguments, see \code{\link[stabs]{stabsel}}. 
#' 
#' @family statistics
#' 
#' @examples
#' p=100; n=100; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#' knockoffs = function(X) create.gaussian(X, mu, Sigma)
#' 
#' # Basic usage with default arguments
#' result = knockoff.filter(X, y, knockoffs=knockoffs,
#'                            statistic=stat.stability_selection)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = stat.stability_selection
#' k_stat = function(X, X_k, y) foo(X, X_k, y, fitfun=stabs::lars.lasso)
#' result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname stat.stability_selection
#' @export
stat.stability_selection <- function(X, X_k, y, fitfun = stabs::glmnet.lasso, ...) {
  if (!requireNamespace('stabs', quietly=T))
    stop('stabs is not installed', call.=F)
  if (!is.vector(y)) {
    stop('Knockoff statistic stat.stability_selection requires the input y to be a vector')
  }
  Z = stability_selection_importance(cbind(X, X_k), y, fitfun=fitfun, ...)
  p = ncol(X)
  orig = 1:p
  abs(Z[orig]) - abs(Z[orig+p])
}

#' Stability selection
#' 
#' Perform variable selection with stability selection
#' 
#' @param X matrix of predictors
#' @param y response vector
#' @return vector with jth component the selection probability of variable j
#' 
#' @keywords internal
stability_selection_importance <- function(X, y, ...) {
  X = normc(X)
  
  if (!methods::hasArg(cutoff) ) {
    cutoff = 0.75
  }
  if (!methods::hasArg(PFER) ) {
    PFER = 1
  }
  
  stabFit = stabs::stabsel(X, y, cutoff=cutoff, PFER=PFER, ...)
  rowMeans(unname(stabFit$phat))
}