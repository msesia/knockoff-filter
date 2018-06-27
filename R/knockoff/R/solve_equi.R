#' Optimization for equi-correlated fixed-X and Gaussian knockoffs
#' 
#' This function solves a very simple optimization problem needed to create fixed-X and 
#' Gaussian SDP knockoffs on the full the covariance matrix. This may be significantly
#' less powerful than \code{\link{create.solve_sdp}}.
#' 
#' @param Sigma positive-definite p-by-p covariance matrix.
#' @return The solution \eqn{s} to the optimization problem defined above.
#' 
#' @details Computes the closed-form solution to the semidefinite programming problem:
#'  \deqn{ \mathrm{maximize}  \; s \quad
#'        \mathrm{subject} \; \mathrm{to:}   \; 0 \leq s \leq 1, \;
#'        2\Sigma - sI \geq 0 }
#' used to generate equi-correlated knockoffs.
#' 
#' The closed form-solution to this problem is \eqn{s = 2\lambda_{\mathrm{min}}(\Sigma) \land 1}.
#' 
#' @family optimization
#' 
#' @export
create.solve_equi <- function(Sigma) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  p = nrow(Sigma)
  tol = 1e-10
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    stop('The covariance matrix is not positive-definite: cannot solve SDP',immediate.=T)
  }
  
  if (p>2) {
    converged=FALSE
    maxitr=10000
    while (!converged) {
      lambda_min = RSpectra::eigs(G, 1, which="SR", opts=list(retvec = FALSE, maxitr=100000, tol=1e-8))$values
      if (length(lambda_min)==1) {
        converged = TRUE
      } else {
        if (maxitr>1e8) {
          warning('In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the 
                covariance matrix. RSpectra::eigs did not converge. Giving up and computing full SVD with built-in R function.',immediate.=T)
          lambda_min = eigen(G, symmetric=T, only.values = T)$values[p]
          converged=TRUE
        } else {
          warning('In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the 
                covariance matrix. RSpectra::eigs did not converge. Trying again with increased number of iterations.',immediate.=T)
          maxitr = maxitr*10
        }
      }
    }
  } else {
    lambda_min = eigen(G, symmetric=T, only.values = T)$values[p]
  }
  
  if (lambda_min<0) {
    stop('In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the 
                covariance matrix. The covariance matrix is not positive-definite.')
  }
  
  s = rep(1, nrow(Sigma)) * min(2*lambda_min, 1)
  
  # Compensate for numerical errors (feasibility)
  psd = 0;
  s_eps = 1e-8;
  while (psd==0) {
    psd = is_posdef(2*G-diag(s*(1-s_eps),length(s)))
    if (!psd) {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)
  
  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}