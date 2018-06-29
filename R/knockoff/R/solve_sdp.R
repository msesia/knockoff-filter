#' Optimization for fixed-X and Gaussian knockoffs
#'
#' This function solves the optimization problem needed to create fixed-X and Gaussian SDP knockoffs
#' on the full covariance matrix. This will be more powerful than \code{\link{create.solve_asdp}},
#' but more computationally expensive.
#' 
#' @param Sigma positive-definite p-by-p covariance matrix.
#' @param maxit maximum number of iterations for the solver (default: 1000).
#' @param gaptol tolerance for duality gap as a fraction of the value of the objective functions (default: 1e-6).
#' @param verbose whether to display progress (default: FALSE).
#' @return The solution \eqn{s} to the semidefinite programming problem defined above.
#'
#' @details
#' Solves the semidefinite programming problem:
#'
#'   \deqn{ \mathrm{maximize}      \; \mathrm{sum}(s) \quad
#'           \mathrm{subject} \; \mathrm{to}    0 \leq s \leq 1, \;
#'                                  2\Sigma - \mathrm{diag}(s) \geq 0}
#' 
#' This problem is solved using the interior-point method implemented in \code{\link[Rdsdp]{dsdp}}.
#'
#' If the matrix Sigma supplied by the user is a non-scaled covariance matrix 
#' (i.e. its diagonal entries are not all equal to 1), then the appropriate scaling is applied before
#' solving the SDP defined above. The result is then scaled back before being returned, as to match 
#' the original scaling of the covariance matrix supplied by the user.
#' 
#' @family optimization
#' 
#' @export
create.solve_sdp <- function(Sigma, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]
  
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }
  
  # Convert problem for SCS
  
  # Linear constraints
  Cl1 = rep(0,p)
  Al1 = -Matrix::Diagonal(p)
  Cl2 = rep(1,p)
  Al2 = Matrix::Diagonal(p)
  
  # Positive-definite cone
  d_As = c(diag(p))
  As = Matrix::Diagonal(length(d_As), x=d_As)
  As = As[which(Matrix::rowSums(As) > 0),] 
  Cs = c(2*G)
  
  # Assemble constraints and cones
  A = cbind(Al1,Al2,As)
  C = matrix(c(Cl1,Cl2,Cs),1)
  K=NULL
  K$s=p
  K$l=2*p
  
  # Objective
  b = rep(1,p)
  
  # Solve SDP with Rdsdp
  OPTIONS=NULL
  OPTIONS$gaptol=gaptol
  OPTIONS$maxit=maxit
  OPTIONS$logsummary=0
  OPTIONS$outputstats=0
  OPTIONS$print=0
  if(verbose) cat("Solving SDP ... ")
  sol = Rdsdp::dsdp(A,b,C,K,OPTIONS)
  if(verbose) cat("done. \n")
  
  # Check whether the solution is feasible
  if( ! identical(sol$STATS$stype,"PDFeasible")) {
    warning('The SDP solver returned a non-feasible solution. Knockoffs may lose power.')
  }
  
  # Clip solution to correct numerical errors (domain)
  s = sol$y
  s[s<0]=0
  s[s>1]=1
  
  # Compensate for numerical errors (feasibility)
  if(verbose) cat("Verifying that the solution is correct ... ")
  psd = 0
  s_eps = 1e-8
  while ((psd==0) & (s_eps<=0.1)) {
    if (is_posdef(2*G-diag(s*(1-s_eps),length(s)),tol=1e-9)) {
      psd  = 1
    }
    else {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)
  s[s<0]=0
  if(verbose) cat("done. \n")
  
  # Verify that the solution is correct
  if (all(s==0)) {
    warning('In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
  }
  
  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}

#' Vectorize a matrix into the SCS format
#'  
#' @rdname vectorize_matrix
#' @keywords internal
create.vectorize_matrix = function(M) {
  # Scale the off-diagonal entries by sqrt(2)
  vectorized_matrix = M
  vectorized_matrix[lower.tri(M,diag=FALSE)] = M[lower.tri(M,diag=FALSE)] * sqrt(2)
  # Stack the lower triangular elements column-wise
  vectorized_matrix = vectorized_matrix[lower.tri(vectorized_matrix,diag=TRUE)]
}