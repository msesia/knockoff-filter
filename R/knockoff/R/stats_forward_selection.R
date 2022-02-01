#' Importance statistics based on forward selection
#' 
#' Computes the statistic
#'   \deqn{W_j = \max(Z_j, Z_{j+p}) \cdot \mathrm{sgn}(Z_j - Z_{j+p}),}
#' where \eqn{Z_1,\dots,Z_{2p}} give the reverse order in which the 2p
#' variables (the originals and the knockoffs) enter the forward selection 
#' model.
#' See the Details for information about forward selection.
#' 
#' In \emph{forward selection}, the variables are chosen iteratively to maximize
#' the inner product with the residual from the previous step. The initial
#' residual is always \code{y}. In standard forward selection
#' (\code{stat.forward_selection}), the next residual is the remainder after
#' regressing on the selected variable; when orthogonal matching pursuit
#' is used, the next residual is the remainder
#' after regressing on \emph{all} the previously selected variables.
#' 
#' @param X    n-by-p matrix of original variables.
#' @param X_k  n-by-p matrix of knockoff variables.
#' @param y    numeric vector of length n, containing the response variables.
#' @param omp  whether to use orthogonal matching pursuit (default: F).
#' @return A vector of statistics \eqn{W} of length p.
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
#'                            statistic=stat.forward_selection)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = stat.forward_selection
#' k_stat = function(X, X_k, y) foo(X, X_k, y, omp=TRUE)
#' result = knockoff.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname stat.forward_selection
#' @export
stat.forward_selection <- function(X, X_k, y, omp=F) {
  if( is.numeric(y) ){
    y = as.vector(y)
  } else {
    stop('Knockoff statistic stat.forward_selection requires the input y to be a numeric vector')
  }
  p = ncol(X)
  X = scale(X)
  X_k = scale(X_k)
  
  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)
  
  # Compute statistics
  path = fs(cbind(X.swap, Xk.swap), y, omp)
  Z = 2*p + 1 - order(path) # Are we recycling here?
  orig = 1:p
  W = pmax(Z[orig], Z[orig+p]) * sign(Z[orig] - Z[orig+p])
  
  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
}

#' Forward selection
#' 
#' Perform forward variable selection with or without OMP
#' 
#' @param X matrix of predictors
#' @param y response vector
#' @param omp whether to use orthogonal matching pursuit (OMP)
#' @return vector with jth component the variable added at step j
#' 
#' @keywords internal
fs <- function(X, y, omp=FALSE) {
  n = nrow(X); p = ncol(X)
  stopifnot(n == length(y))
  path = rep.int(0, p)
  in_model = rep(FALSE, p)
  residual = y
  if (omp) Q = matrix(0, n, p)
  
  for (step in 1:p) {
    # Find the best variable to add among the remaining variables.
    available_vars = which(!in_model)
    products = apply(X[,!in_model,drop=F], 2,
                     function(x) abs(sum(x * residual)))
    best_var = available_vars[which.max(products)][1]
    path[step] = best_var
    in_model[best_var] = TRUE
    
    # Update the residual.
    x = X[,best_var]
    if (step == p) break
    if (omp) {
      for (j in seq(1, length.out=step-1))
        x = x - Q[,j]%*%x * Q[,j]
      q = x / sqrt(sum(x^2))
      Q[,step] = q
      residual = residual - q%*%y * q
    } 
    else {
      residual = residual - as.vector(x %*% residual) * as.vector(x)
    }
  }
  return(path)
}
