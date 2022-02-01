#' Importance statistics based on random forests
#' 
#' Computes the difference statistic
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the random forest feature importances
#' of the jth variable and its knockoff, respectively.
#' 
#' @param X n-by-p matrix of original variables.
#' @param X_k n-by-p matrix of knockoff variables.
#' @param y vector of length n, containing the response variables. If a factor, classification is assumed, 
#' otherwise regression is assumed.
#' @param ... additional arguments specific to \code{ranger} (see Details).
#' @return A vector of statistics \eqn{W} of length p.
#' 
#' @details This function uses the \code{ranger} package to compute variable 
#' importance measures. The importance of a variable is measured as the total decrease
#' in node impurities from splitting on that variable, averaged over all trees. 
#' For regression, the node impurity is measured by residual sum of squares.
#' For classification, it is measured by the Gini index.
#' 
#' For a complete list of the available additional arguments, see \code{\link[ranger]{ranger}}. 
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
#' y = X %*% beta + rnorm(n)
#' knockoffs = function(X) create.gaussian(X, mu, Sigma)
#' 
#' # Basic usage with default arguments
#' result = knockoff.filter(X, y, knockoffs=knockoffs, 
#'                            statistic=stat.random_forest)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = stat.random_forest
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nodesize=5)
#' result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname stat.random_forest
#' @export
stat.random_forest <- function(X, X_k, y, ...) {
  if (!requireNamespace('ranger', quietly=T))
    stop('ranger is not installed', call.=F)
  
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)
  
  # Compute statistics
  Z = random_forest_importance(cbind(X.swap, Xk.swap), y) 
  p = ncol(X)
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig+p])
  
  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
}

#' @keywords internal
random_forest_importance <- function(X, y, ...) {
  df = data.frame(y=y, X=X)
  rfFit = ranger::ranger(y~., data=df, importance="impurity", write.forest=F, ...)
  as.vector(rfFit$variable.importance)
}