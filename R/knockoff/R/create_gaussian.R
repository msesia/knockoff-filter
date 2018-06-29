#' Model-X Gaussian knockoffs
#' 
#' Samples multivariate Gaussian model-X knockoff variables.
#' 
#' @param X n-by-p matrix of original variables.
#' @param mu vector of length p, indicating the mean parameter of the Gaussian model for \eqn{X}.
#' @param Sigma p-by-p covariance matrix for the Gaussian model of \eqn{X}.
#' @param method either "equi", "sdp" or "asdp" (default: "asdp").
#' This determines the method that will be used to minimize the correlation between the original variables and the knockoffs.
#' @param diag_s vector of length p, containing the pre-computed covariances between the original 
#' variables and the knockoffs. This will be computed according to \code{method}, if not supplied. 
#' @return A n-by-p matrix of knockoff variables.
#' 
#' @family create
#' 
#' @references 
#'   Candes et al., Panning for Gold: Model-free Knockoffs for High-dimensional Controlled Variable Selection,
#'   arXiv:1610.02351 (2016).
#'   \href{https://web.stanford.edu/group/candes/knockoffs/index.html}{https://web.stanford.edu/group/candes/knockoffs/index.html}
#' 
#' @examples
#' p=200; n=100; k=15
#' rho = 0.4
#' mu = rep(0,p); Sigma = toeplitz(rho^(0:(p-1)))
#' X = matrix(rnorm(n*p),n) %*% chol(Sigma)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#' 
#' # Basic usage with default arguments
#' knockoffs = function(X) create.gaussian(X, mu, Sigma)
#' result = knockoff.filter(X, y, knockoffs=knockoffs)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' knockoffs = function(X) create.gaussian(X, mu, Sigma, method='equi')
#' result = knockoff.filter(X, y, knockoffs=knockoffs)
#' print(result$selected)
#' 
#' @export
create.gaussian <- function(X, mu, Sigma, method=c("asdp","sdp","equi"), diag_s=NULL) {
  method = match.arg(method)
  
  # Do not use ASDP unless p>500
  if ((nrow(Sigma)<=500) && method=="asdp") {
    method="sdp"
  }
  
  if (is.null(diag_s)) {
    diag_s = diag(switch(match.arg(method),
                    'equi' = create.solve_equi(Sigma),
                    'sdp'  = create.solve_sdp(Sigma),
                    'asdp' = create.solve_asdp(Sigma)))
  }
  if (is.null(dim(diag_s))) {
    diag_s = diag(diag_s,length(diag_s))
  }
  
  # If diag_s is zero, we can only generate trivial knockoffs.
  if(all(diag_s==0)) {
    warning("The conditional knockoff covariance matrix is not positive definite. Knockoffs will have no power.")
    return(X)
  }
  
  SigmaInv_s = solve(Sigma,diag_s)
  mu_k = X - sweep(X,2,mu,"-") %*% SigmaInv_s
  Sigma_k = 2*diag_s - diag_s %*% SigmaInv_s
  X_k = mu_k + matrix(rnorm(ncol(X)*nrow(X)),nrow(X)) %*% chol(Sigma_k)
}